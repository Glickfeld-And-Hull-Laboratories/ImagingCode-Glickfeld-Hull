"""Central variable store: MWEL defaults + XML overrides."""

from dataclasses import dataclass
from mwel_parser import parse_mwel, MwelVariable, MwelSelection
from variable_set_loader import load_variable_set


@dataclass
class Variable:
    """Runtime variable with current value and metadata."""
    name: str
    var_type: str
    default_value: object
    current_value: object
    groups: str = ""
    enclosing_group: str = ""
    persistent: bool = False
    overridden: bool = False  # True if XML override was applied

    @property
    def is_modified(self) -> bool:
        return self.current_value != self.default_value


class ExperimentModel:
    """Unified variable store for the MWorks trial simulator.

    Loads defaults from .mwel, applies overrides from .xml variable-set files.
    Protocol engine calls get()/set() during trial execution.
    """

    def __init__(self, mwel_path: str, xml_override_path: str | None = None):
        self.mwel_path = mwel_path
        self.xml_override_path = xml_override_path

        # Parse MWEL
        mwel_vars, mwel_sels = parse_mwel(mwel_path)
        self.mwel_selections: dict[str, MwelSelection] = mwel_sels

        # Build variable store from MWEL defaults
        self.variables: dict[str, Variable] = {}
        self.groups: dict[str, list[str]] = {}  # group name -> sorted var names

        for name, mv in mwel_vars.items():
            var = Variable(
                name=name,
                var_type=mv.var_type,
                default_value=mv.default_value,
                current_value=mv.default_value,
                groups=mv.groups,
                enclosing_group=mv.enclosing_group,
                persistent=mv.persistent,
            )
            self.variables[name] = var

            # Index by groups attribute
            if mv.groups:
                self.groups.setdefault(mv.groups, []).append(name)

        # Sort group lists for stable display
        for g in self.groups:
            self.groups[g].sort()

        # Apply XML overrides if provided
        if xml_override_path:
            self.load_overrides(xml_override_path)

    def get(self, name: str, default=None):
        """Get current value of a variable."""
        var = self.variables.get(name)
        if var is None:
            return default
        return var.current_value

    def get_float(self, name: str, default: float = 0.0) -> float:
        """Get current value as float."""
        val = self.get(name)
        if val is None:
            return default
        try:
            return float(val)
        except (ValueError, TypeError):
            return default

    def get_int(self, name: str, default: int = 0) -> int:
        """Get current value as int."""
        val = self.get(name)
        if val is None:
            return default
        try:
            return int(float(val))
        except (ValueError, TypeError):
            return default

    def get_bool(self, name: str, default: bool = False) -> bool:
        """Get current value as bool."""
        val = self.get(name)
        if val is None:
            return default
        if isinstance(val, bool):
            return val
        if isinstance(val, (int, float)):
            return bool(val)
        if isinstance(val, str):
            return val.lower() not in ("0", "false", "no", "")
        return bool(val)

    def set(self, name: str, value):
        """Set current value of a variable. Creates it if it doesn't exist."""
        if name in self.variables:
            self.variables[name].current_value = value
        else:
            # Create new variable (for t-vars set by protocol engine)
            self.variables[name] = Variable(
                name=name,
                var_type="float" if isinstance(value, float) else
                         "integer" if isinstance(value, int) else
                         "bool" if isinstance(value, bool) else "string",
                default_value=value,
                current_value=value,
            )

    def reset(self, name: str):
        """Reset a variable to its default value."""
        if name in self.variables:
            self.variables[name].current_value = self.variables[name].default_value

    def reset_to_defaults(self):
        """Reset all variables to MWEL defaults (removes all overrides)."""
        for var in self.variables.values():
            var.current_value = var.default_value
            var.overridden = False

    def load_overrides(self, xml_path: str):
        """Apply variable overrides from an XML variable-set file."""
        self.xml_override_path = xml_path
        overrides = load_variable_set(xml_path)

        for name, (type_str, value) in overrides.items():
            if name in self.variables:
                self.variables[name].current_value = value
                self.variables[name].overridden = True
            else:
                # Variable not in MWEL — add it anyway
                self.variables[name] = Variable(
                    name=name,
                    var_type=type_str,
                    default_value=value,
                    current_value=value,
                    overridden=True,
                )

    def get_modified_count(self) -> int:
        """Count variables that differ from their defaults."""
        return sum(1 for v in self.variables.values() if v.is_modified)

    def get_overridden_count(self) -> int:
        """Count variables that were overridden by XML."""
        return sum(1 for v in self.variables.values() if v.overridden)

    @property
    def enclosing_groups(self) -> dict[str, list[str]]:
        """Group vars by their enclosing group { } block."""
        result: dict[str, list[str]] = {}
        for name, var in self.variables.items():
            eg = var.enclosing_group or "Ungrouped"
            result.setdefault(eg, []).append(name)
        for g in result:
            result[g].sort()
        return result
