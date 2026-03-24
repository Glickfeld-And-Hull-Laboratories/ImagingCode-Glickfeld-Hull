"""Parse variable and selection declarations from .mwel files."""

import re
from dataclasses import dataclass, field


@dataclass
class MwelVariable:
    """A variable declaration parsed from .mwel."""
    name: str
    var_type: str = "float"        # bool, float, integer, string
    default_value: object = 0
    persistent: bool = False
    groups: str = ""
    enclosing_group: str = ""      # group 'Name' { } block this var is inside


@dataclass
class MwelSelection:
    """A selection variable declaration parsed from .mwel."""
    name: str
    values: list = field(default_factory=lambda: list(range(80)))
    selection_type: str = "random_without_replacement"
    nsamples: int = 80
    sampling_method: str = "samples"


# Regex patterns for variable declarations
# Matches: var name = (type)(default) ...
_VAR_TYPED_DEFAULT = re.compile(
    r'var\s+(\w+)\s*=\s*\((\w+)\)\s*\(([^)]*)\)'
)
# Matches: var name = default (no type cast) — numeric or string literal
_VAR_BARE_DEFAULT = re.compile(
    r'var\s+(\w+)\s*=\s*(?!\()([^\s(]+)'
)
# Matches: var name = [list,items]
_VAR_LIST_DEFAULT = re.compile(
    r'var\s+(\w+)\s*=\s*\[([^\]]*)\]'
)
# Matches: groups = 'Group Name' or groups = GroupName
_GROUPS_ATTR = re.compile(
    r"groups\s*=\s*(?:'([^']*)'|(\w+))"
)
# Matches: persistant = 1
_PERSISTENT_ATTR = re.compile(
    r'persistant\s*=\s*1'
)
# Matches: group 'Name' { or group Name {
_GROUP_BLOCK = re.compile(
    r"group\s+(?:'([^']*)'|(\w+))\s*\{")
# Matches: selection svName (
_SELECTION_DECL = re.compile(
    r'selection\s+(\w+)\s*\(')
# Matches: values = 0,1,2,...
_SEL_VALUES = re.compile(
    r'values\s*=\s*([\d,\s]+)')
# Matches: nsamples = N
_SEL_NSAMPLES = re.compile(
    r'nsamples\s*=\s*(\d+)')
# Matches: selection = random_without_replacement
_SEL_TYPE = re.compile(
    r'selection\s*=\s*(\w+)')


def _coerce_value(val_str: str, type_str: str = ""):
    """Convert a string value to Python type based on declared type."""
    val_str = val_str.strip().strip("'\"")

    if type_str == "bool":
        return val_str not in ("0", "false", "False", "NO")
    if type_str == "float":
        try:
            return float(val_str)
        except ValueError:
            return 0.0
    if type_str == "integer":
        try:
            return int(float(val_str))
        except ValueError:
            return 0
    if type_str == "string":
        return val_str

    # No type declared — guess from value
    if val_str.lower() in ("true", "false", "yes", "no"):
        return val_str.lower() in ("true", "yes")
    try:
        if '.' in val_str:
            return float(val_str)
        return int(val_str)
    except ValueError:
        return val_str


def _infer_type(val_str: str) -> str:
    """Infer MWEL type from a bare value string."""
    val_str = val_str.strip()
    if val_str.lower() in ("true", "false", "yes", "no"):
        return "bool"
    try:
        if '.' in val_str:
            float(val_str)
            return "float"
        int(val_str)
        return "integer"
    except ValueError:
        return "string"


def parse_mwel(filepath: str) -> tuple[dict[str, MwelVariable], dict[str, MwelSelection]]:
    """Parse a .mwel file and return (variables, selections).

    Returns:
        variables: dict mapping variable name to MwelVariable
        selections: dict mapping selection name to MwelSelection
    """
    with open(filepath, 'r') as f:
        text = f.read()

    lines = text.split('\n')
    variables: dict[str, MwelVariable] = {}
    selections: dict[str, MwelSelection] = {}

    # Track enclosing group blocks via brace depth
    group_stack: list[str] = []
    brace_depth = 0
    group_start_depth: list[int] = []

    # Track attribute context for multi-line var declarations
    current_var: MwelVariable | None = None
    in_var_attrs = False
    paren_depth = 0

    # Track selection context
    current_sel: MwelSelection | None = None
    in_sel_attrs = False
    sel_paren_depth = 0

    for line in lines:
        stripped = line.strip()

        # Skip comments
        if stripped.startswith('//'):
            # Count braces in comments? No, skip entirely
            pass
        else:
            # Check for group block opening
            gm = _GROUP_BLOCK.search(stripped)
            if gm:
                group_name = gm.group(1) or gm.group(2)
                group_stack.append(group_name)
                # Count braces on this line
                open_b = stripped.count('{')
                close_b = stripped.count('}')
                group_start_depth.append(brace_depth + open_b)
                brace_depth += open_b - close_b
                # Check for selection on this line too
                sm = _SELECTION_DECL.search(stripped)
                if sm:
                    current_sel = MwelSelection(name=sm.group(1))
                    in_sel_attrs = True
                    sel_paren_depth = stripped.count('(') - stripped.count(')')
                continue

            # Track brace depth for group blocks (but not inside var/selection parens)
            if not in_var_attrs and not in_sel_attrs:
                open_b = stripped.count('{')
                close_b = stripped.count('}')
                brace_depth += open_b - close_b
                # Pop group stack when depth returns
                while group_start_depth and brace_depth < group_start_depth[-1]:
                    group_start_depth.pop()
                    if group_stack:
                        group_stack.pop()

            # Check for selection declaration
            sm = _SELECTION_DECL.search(stripped)
            if sm and not in_var_attrs:
                current_sel = MwelSelection(name=sm.group(1))
                in_sel_attrs = True
                sel_paren_depth = stripped.count('(') - stripped.count(')')
                # Parse inline attributes
                _parse_sel_attrs(stripped, current_sel)
                if sel_paren_depth <= 0:
                    selections[current_sel.name] = current_sel
                    current_sel = None
                    in_sel_attrs = False
                continue

            # Continue parsing selection attributes
            if in_sel_attrs and current_sel:
                sel_paren_depth += stripped.count('(') - stripped.count(')')
                _parse_sel_attrs(stripped, current_sel)
                if sel_paren_depth <= 0:
                    selections[current_sel.name] = current_sel
                    current_sel = None
                    in_sel_attrs = False
                continue

            # Check for variable declarations
            # List default: var name = [items]
            lm = _VAR_LIST_DEFAULT.search(stripped)
            if lm:
                name = lm.group(1)
                items_str = lm.group(2)
                items = [_coerce_value(x.strip()) for x in items_str.split(',') if x.strip()]
                var = MwelVariable(
                    name=name,
                    var_type="list",
                    default_value=items,
                    enclosing_group=group_stack[-1] if group_stack else ""
                )
                current_var = var
                in_var_attrs = True
                paren_depth = stripped.count('(') - stripped.count(')')
                _parse_var_attrs(stripped, current_var, lm.end())
                if paren_depth <= 0 and ')' not in stripped[lm.end():]:
                    variables[var.name] = var
                    current_var = None
                    in_var_attrs = False
                continue

            # Typed default: var name = (type)(default)
            tm = _VAR_TYPED_DEFAULT.search(stripped)
            if tm:
                name = tm.group(1)
                vtype = tm.group(2)
                default = _coerce_value(tm.group(3), vtype)
                var = MwelVariable(
                    name=name,
                    var_type=vtype,
                    default_value=default,
                    enclosing_group=group_stack[-1] if group_stack else ""
                )
                current_var = var
                in_var_attrs = True
                paren_depth = stripped.count('(') - stripped.count(')')
                _parse_var_attrs(stripped, current_var, tm.end())
                if paren_depth <= 0:
                    variables[var.name] = var
                    current_var = None
                    in_var_attrs = False
                continue

            # Bare default: var name = value
            bm = _VAR_BARE_DEFAULT.search(stripped)
            if bm:
                name = bm.group(1)
                val_str = bm.group(2).strip()
                # Skip if val_str starts with ( — it's a typed default we missed
                if val_str.startswith('('):
                    continue
                vtype = _infer_type(val_str)
                default = _coerce_value(val_str, vtype)
                var = MwelVariable(
                    name=name,
                    var_type=vtype,
                    default_value=default,
                    enclosing_group=group_stack[-1] if group_stack else ""
                )
                current_var = var
                in_var_attrs = True
                paren_depth = stripped.count('(') - stripped.count(')')
                _parse_var_attrs(stripped, current_var, bm.end())
                if paren_depth <= 0:
                    variables[var.name] = var
                    current_var = None
                    in_var_attrs = False
                continue

            # Continue parsing variable attributes
            if in_var_attrs and current_var:
                paren_depth += stripped.count('(') - stripped.count(')')
                _parse_var_attrs(stripped, current_var, 0)
                if paren_depth <= 0:
                    variables[current_var.name] = current_var
                    current_var = None
                    in_var_attrs = False
                continue

    # Finalize any remaining var/selection
    if current_var:
        variables[current_var.name] = current_var
    if current_sel:
        selections[current_sel.name] = current_sel

    return variables, selections


def _parse_var_attrs(line: str, var: MwelVariable, start_pos: int = 0):
    """Extract groups and persistant attributes from a line."""
    region = line[start_pos:]
    gm = _GROUPS_ATTR.search(region)
    if gm:
        var.groups = gm.group(1) or gm.group(2) or ""
    if _PERSISTENT_ATTR.search(region):
        var.persistent = True


def _parse_sel_attrs(line: str, sel: MwelSelection):
    """Extract selection attributes from a line."""
    vm = _SEL_VALUES.search(line)
    if vm:
        vals_str = vm.group(1)
        sel.values = [int(x.strip()) for x in vals_str.split(',') if x.strip()]
    nm = _SEL_NSAMPLES.search(line)
    if nm:
        sel.nsamples = int(nm.group(1))
    # selection = type (but not the outer 'selection svName')
    for m in _SEL_TYPE.finditer(line):
        val = m.group(1)
        if val != sel.name and val in ("random_without_replacement", "sequential", "random_with_replacement"):
            sel.selection_type = val
