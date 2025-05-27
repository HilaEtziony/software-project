from enum import Enum, auto
from pathlib import Path
import re
import subprocess
import sys

LINE_CONFIG_REGEX = re.compile(
    r"(?P<idx>\d+)\. k=(?P<k>\d+), max_iter\s+=\s+(?:not provided|(?P<max_iter>\d+))"
)
DEFAULT_MAX_ITER = "400"
TEST_PYTHON = True
TEST_C = True

class TestType(Enum):
    PYTHON = auto()
    C = auto()


def main():
    tests_dir = Path(sys.argv[1])

    readme_path = tests_dir / Path("test_readme.txt")
    configs = []
    with readme_path.open() as f:
        for line in f:
            match = LINE_CONFIG_REGEX.match(line)
            if match:
                config = match.groupdict()
                if config["max_iter"] is None:
                    config["max_iter"] = DEFAULT_MAX_ITER
                configs.append(config)

    for config in configs:
        print("Test", config["idx"])
        if TEST_PYTHON:
            print("testing python")
            run_test(tests_dir, config, TestType.PYTHON)

        if TEST_C:
            print("testing C")
            run_test(tests_dir, config, TestType.PYTHON)
        
        if TEST_PYTHON or TEST_C:
            print()


def run_test(tests_dir: Path, config: dict[str, str], test_type: TestType):
    success = True

    input_path = tests_dir / Path(f"input_{config['idx']}.txt")
    output_path = tests_dir / Path(f"output_{config['idx']}.txt")
    with input_path.open() as input_file:
        match test_type:
            case TestType.PYTHON:
                args =["python3", "kmeans.py", config["k"], config["max_iter"]]
            case TestType.C:
                args =["kmeans", config["k"], config["max_iter"]]
    
        result = subprocess.run(
                args,
                stdin=input_file,
                capture_output=True,
                text=True,
            )

    if result.returncode:
        print(f"process returned with code {result.returncode}")
        success = False

    if result.stderr:
        print("process had non-empty stderr:")
        print(result.stderr)
        success = False

        # Compare outputs
    if output_path.read_text().rstrip() != result.stdout.rstrip():
        print("process had mismatching output")
        success = False
        
    if success:
        print("\033[32msuccess\033[0m")
    else:
        print("\033[31mfailure\033[0m")


if __name__ == "__main__":
    main()
