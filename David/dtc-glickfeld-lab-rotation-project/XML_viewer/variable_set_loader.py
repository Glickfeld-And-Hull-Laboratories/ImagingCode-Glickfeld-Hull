"""Parse variable_assignments from MWorks XML variable-set files."""

import xml.etree.ElementTree as ET


def _coerce_xml_value(type_str: str, value_str: str):
    """Convert XML value string to Python type."""
    t = type_str.lower()
    if t == "boolean":
        return value_str.lower() not in ("false", "0", "no")
    if t == "integer":
        return int(float(value_str))
    if t == "float":
        return float(value_str)
    if t == "string":
        return value_str
    # Fallback: try numeric
    try:
        if '.' in value_str:
            return float(value_str)
        return int(value_str)
    except ValueError:
        return value_str


def load_variable_set(filepath: str) -> dict[str, tuple[str, object]]:
    """Load an XML variable-set file.

    Returns:
        dict mapping variable name to (type_string, coerced_value).
        For list variables, type_string is "list".
    """
    tree = ET.parse(filepath)
    root = tree.getroot()

    # Find <variable_assignments> element
    va_elem = root.find('variable_assignments')
    if va_elem is None:
        va_elem = root  # try root directly

    result: dict[str, tuple[str, object]] = {}

    for elem in va_elem.findall('variable_assignment'):
        name = elem.get('variable')
        if not name:
            continue

        type_str = elem.get('type', '')
        value_str = elem.get('value', '')

        # Check for list_data child (e.g., LED_sequence)
        list_data = elem.find('list_data')
        if list_data is not None:
            items = []
            for list_elem in list_data.findall('list_element'):
                val_elem = list_elem.find('value')
                if val_elem is not None:
                    item_type = val_elem.get('type', 'integer')
                    item_val = val_elem.get('value', '0')
                    items.append(_coerce_xml_value(item_type, item_val))
            result[name] = ("list", items)
        elif type_str and value_str is not None:
            result[name] = (type_str, _coerce_xml_value(type_str, value_str))
        elif value_str is not None:
            # No type specified, try to infer
            result[name] = ("string", value_str)

    return result
