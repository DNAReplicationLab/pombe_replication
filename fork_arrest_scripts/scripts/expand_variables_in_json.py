import json
import sys

# Script written using ChatGPT


def expand_variables_in_json(input_json):
    """    Perform replacements in the input json file according to the other keys and values available
    Args:
        input_json: python dictionary representing the input json file

    Returns:
        python dictionary representing the processed json file with replacements made

    Example:
        >>> test_json = {
        ...   "key1": "blah_${var1}_blah",
        ...   "key2": "${key1}_ipsum",
        ...   "key3": "blah",
        ...   "key4": 123,
        ...   "var1": "value1",
        ...   "var2": "value2"
        ... }
        >>> expand_variables_in_json(test_json)
        {'key1': 'blah_value1_blah', 'key2': 'blah_value1_blah_ipsum', 'key3': 'blah', 'key4': 123, 'var1': 'value1', 'var2': 'value2'}
    """
    # Construct a dictionary using all key-value pairs where value is a string
    string_values = {key: value for key, value in input_json.items() if isinstance(value, str)}

    # Process the input JSON and perform replacements.
    # NOTE: We are doing it in a lazy way for multiple levels of replacements (i.e. we check for replacements 10 times)
    #       and by brute force checking i.e. for every key-value pair we check every other key-value pair
    #       It is possible to improve this logic, but we don't care for now.
    for _ in range(10):
        for key, value in input_json.items():
            if isinstance(value, str):
                for string_key, string_value in string_values.items():
                    if string_key != key:
                        value = value.replace("${" + string_key + "}", string_value)
                input_json[key] = value

    # Perform a final check for the ${ construction and raise an exception if found
    for key, value in input_json.items():
        if isinstance(value, str):
            if "${" in value:
                raise Exception("Error: Could not perform all replacements in the JSON file. "
                                "Please check the following key-value pair: " + key + ": " + value)

    return input_json


def main_expand_vars_in_json():
    # Read JSON from stdin
    input_data = json.load(sys.stdin)

    # Process JSON
    processed_data = expand_variables_in_json(input_data)

    # Output processed JSON to stdout
    json.dump(processed_data, sys.stdout, indent=4)


if __name__ == "__main__":

    if sys.stdin.isatty():
        print("Usage: cat input.json | python expand_variables_in_json.py > output.json")
        print("Goal: Any key-value pair in the json file where the value is a string and contains ${...} will be"
              " replaced with the value of the key in the json file")
        print("Example input json file:")
        print("{")
        print("    \"key1\": \"blah_${var1}_blah\",")
        print("    \"key2\": \"lorem_${var2}_ipsum\",")
        print("    \"var1\": \"value1\",")
        print("    \"var2\": \"value2\"")
        print("}")
        print("Output:")
        print("{")
        print("    \"key1\": \"blah_value1_blah\",")
        print("    \"key2\": \"lorem_value2_ipsum\",")
        print("    \"var1\": \"value1\",")
        print("    \"var2\": \"value2\"")
        print("}")
        sys.exit(1)

    main_expand_vars_in_json()
