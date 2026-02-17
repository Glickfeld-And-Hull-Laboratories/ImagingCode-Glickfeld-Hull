from dataclasses import dataclass, field
from collections import OrderedDict
from typing import Any, Optional
import xml.etree.ElementTree as ET
import copy
import os

@dataclass
class MWVariable:
    tag: str                    # variable name (e.g. "stimOneGratingContrast")
    var_type: str               # "float", "integer", "boolean", "string", "list", "selection"
    default_value: str          # raw string from XML
    current_value: Any          # coerced Python value
    folder: str                 # parent folder tag (e.g. "Stimuli", "Behavioral Control")
    groups: str = ""            # groups attribute from XML (e.g. "Stimulus One Parameters")
    persistant: bool = False    # whether this is a persistent variable
    is_modified: bool = False   # True if current_value differs from default
    xml_element: Optional[Any] = field(default=None, repr=False)  # ET.Element reference

class ExperimentModel:
    def __init__(self, xml_path: str):
        self.xml_path = xml_path
        self.tree = None
        self.root = None
        self.variables: OrderedDict[str, MWVariable] = OrderedDict()
        self.groups: OrderedDict[str, list] = OrderedDict()       # group_name -> [tag, ...]
        self.folder_vars: OrderedDict[str, list] = OrderedDict()  # folder_name -> [tag, ...]
        self._load()

    def _coerce_value(self, raw: str, var_type: str) -> Any:
        """Coerce raw string to proper Python type."""
        if not raw:
            raw = "0"
        vt = var_type.lower()
        if vt == "float":
            try:
                return float(raw)
            except ValueError:
                return 0.0
        elif vt == "integer":
            try:
                return int(float(raw))
            except ValueError:
                return 0
        elif vt == "boolean":
            return raw not in ('0', 'false', '', 'NO')
        else:  # string, list, selection
            return raw

    def _load(self):
        self.tree = ET.parse(self.xml_path)
        self.root = self.tree.getroot()

        # Find all <variable> elements inside <variables> section
        variables_section = self.root.find('variables')
        if variables_section is None:
            return

        for folder_elem in variables_section:
            if folder_elem.tag != 'folder':
                continue
            folder_name = folder_elem.get('tag', 'Unknown')
            if folder_name not in self.folder_vars:
                self.folder_vars[folder_name] = []

            for var_elem in folder_elem:
                if var_elem.tag != 'variable':
                    continue
                tag = var_elem.get('tag', '')
                if not tag:
                    continue

                var_type = var_elem.get('type', 'string').lower()
                default_val = var_elem.get('default_value', '0')
                groups = var_elem.get('groups', '')
                persistant = var_elem.get('persistant', '0') not in ('0', '', None)

                current = self._coerce_value(default_val, var_type)

                mv = MWVariable(
                    tag=tag,
                    var_type=var_type,
                    default_value=default_val,
                    current_value=current,
                    folder=folder_name,
                    groups=groups,
                    persistant=persistant,
                    is_modified=False,
                    xml_element=var_elem,
                )
                self.variables[tag] = mv
                self.folder_vars[folder_name].append(tag)

                if groups:
                    if groups not in self.groups:
                        self.groups[groups] = []
                    self.groups[groups].append(tag)

    def get(self, tag: str) -> Any:
        """Get the current value of a variable."""
        if tag in self.variables:
            return self.variables[tag].current_value
        return None

    def set(self, tag: str, value: Any):
        """Set a variable's value and mark it as modified."""
        if tag not in self.variables:
            return
        mv = self.variables[tag]
        mv.current_value = value
        # Check if modified from default
        default_coerced = self._coerce_value(mv.default_value, mv.var_type)
        mv.is_modified = (mv.current_value != default_coerced)

    def save(self):
        """Write modified default_value attributes back to XML."""
        for tag, mv in self.variables.items():
            if mv.is_modified and mv.xml_element is not None:
                # Convert current_value back to string for XML
                if mv.var_type == 'boolean':
                    mv.xml_element.set('default_value', '1' if mv.current_value else '0')
                else:
                    mv.xml_element.set('default_value', str(mv.current_value))
                mv.default_value = mv.xml_element.get('default_value')
                mv.is_modified = False

        # Write back to file, preserving format as much as possible
        self.tree.write(self.xml_path, xml_declaration=True, encoding='unicode')

    def reset(self, tag: str):
        """Reset a variable to its default value."""
        if tag not in self.variables:
            return
        mv = self.variables[tag]
        mv.current_value = self._coerce_value(mv.default_value, mv.var_type)
        mv.is_modified = False

    def reset_all(self):
        """Reset all variables to their default values."""
        for tag in self.variables:
            self.reset(tag)

    def get_modified_count(self) -> int:
        """Return count of modified variables."""
        return sum(1 for mv in self.variables.values() if mv.is_modified)
