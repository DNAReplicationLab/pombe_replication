import json
import os
import sys


def is_valid_path(path):
    """Check if the given path corresponds to a valid file or folder."""
    return os.path.isfile(path) or os.path.isdir(path)


def main(value):
    # Read the comma-separated list from the user
    input_list = value.split(",")

    # Read the JSON data from the piped file
    try:
        json_data = json.load(sys.stdin)
    except json.JSONDecodeError:
        print("Invalid JSON input.")
        return

    # Check if any of the keys in the input_list have values that are valid paths
    for key in input_list:
        value = json_data.get(key.strip())
        if value and isinstance(value, str) and (not is_valid_path(value)):
            print(False)
            return

    print(True)


if __name__ == "__main__":

    usage_string = "Usage: cat <json_file> | python check_if_file_exists_in_json.py <comma-separated list of keys>"

    if sys.stdin.isatty():
        print("No JSON input piped.")
        print(usage_string)
        sys.exit(1)

    if len(sys.argv) != 2:
        print("Invalid number of arguments.")
        print(usage_string)
        sys.exit(1)

    main(sys.argv[1])
