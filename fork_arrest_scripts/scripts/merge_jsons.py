import json
import os


def load_json_file(path):
    # check if file exists or raise error
    if not os.path.exists(path):
        raise FileNotFoundError(f"File '{path}' does not exist.")

    with open(path, 'r') as f:
        return json.load(f)


def write_json_file(path, data):
    with open(path, 'w') as f:
        json.dump(data, f, indent=4)


def merge_jsons(json1, json2):
    output = {}

    # Combine fields from both JSONs
    all_keys = set(json1.keys()) | set(json2.keys())

    for key in all_keys:
        val1 = json1.get(key)
        val2 = json2.get(key)

        # If key only exists in one JSON or values are the same
        if val1 == val2 or val2 is None or val1 is None:
            output[key] = val1 if val1 is not None else val2
        else:
            # If there's a conflict, ask the user which value to keep
            print(f"Conflict for key '{key}':")
            print(f"1) {val1}")
            print(f"2) {val2}")

            choice_prompt = "Choose which value to keep (1/2) or choose 3 to enter new value: "
            choice = input(choice_prompt).strip()

            while choice not in ['1', '2', '3']:
                print("Invalid choice. Please choose 1 or 2 or 3.")
                choice = input(choice_prompt).strip()

            if choice == '3':
                new_val = input("Enter new value: ").strip()
                output[key] = new_val
            else:
                output[key] = val1 if choice == '1' else val2

        print("Using value: ", output[key], " for key: ", key)
        choice = input("Add the prefix FIX: to remind you to fix it? (y/other input = no)").strip()
        if choice == 'y':
            output[key] = 'FIX: ' + str(output[key])

    return output


if __name__ == '__main__':
    import sys

    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <input_json1> <input_json2> <output_json>")
        print("Goal: combine two json files, resolving conflicts by asking the user which value to keep.")
        print("Program written using ChatGPT")
        exit(1)

    json1_path, json2_path, output_path = sys.argv[1], sys.argv[2], sys.argv[3]

    json_1 = load_json_file(json1_path)
    json_2 = load_json_file(json2_path)

    merged = merge_jsons(json_1, json_2)

    write_json_file(output_path, merged)

    print("Merged JSON saved to", output_path)
