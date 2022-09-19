from pathlib import Path
import sys

def bigscape_copy(input_index, output_main):
    input = Path(input_index)
    output = Path(output_main).parent
    output.mkdir(parents=True, exist_ok=True)
    for item in input.parent.glob("*"):
        target = output / item.name
        target.symlink_to(item.resolve(), target_is_directory=True)
    return

if __name__ == "__main__":
    bigscape_copy(sys.argv[1], sys.argv[2])