#!/usr/bin/env python3
import argparse
import io
from enum import Enum, auto
from pathlib import Path
import random
import re
import subprocess
import tempfile
from typing import Optional

LINE_CONFIG_REGEX = re.compile(
    r"(?P<idx>\d+)\. k=(?P<k>\d+), max_iter\s+=\s+(?:not provided|(?P<max_iter>\d+))"
)

DUMMY_INPUT = """\
1,2,3
4,5,6
7,8,9
1,2,3

"""

INVALID_PARAMETERS_VALUES = (
    "-1",
    "0",
    "bug",
    "באג",
    "🐞",
    "虫",
    "65,536",
    "65536",
    "2.5",
)


class TestType(Enum):
    PYTHON = auto()
    C = auto()


def print_green(msg: str):
    print(f"\033[32m{msg}\033[0m")


def print_red(msg: str):
    print(f"\033[31m{msg}\033[0m")


def print_white_on_red(msg: str):
    print(f"\033[97;41m{msg}\033[0m")


def print_white_on_blue(msg: str):
    print(f"\033[97;44m{msg}\033[0m")


def setup_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "tests_dir",
        help="path to the tests directory (from Moodle)",
        metavar="TESTS_DIR",
    )
    language_options = parser.add_mutually_exclusive_group()
    language_options.add_argument(
        "--only-python",
        action="store_false",
        help="only test Python code",
        dest="test_c",
    )
    language_options.add_argument(
        "--only-c", action="store_false", help="only test C code", dest="test_python"
    )

    return parser


def generate_random_input():
    file = tempfile.TemporaryFile()
    buf = io.TextIOWrapper(file)
    dim = random.randint(2, 10)
    N = random.randint(100, 700)

    centroids = tuple(
        tuple(random.uniform(-20, 20) for _ in range(dim)) for _ in range(5)
    )

    for _ in range(N):
        centroid = random.choice(centroids)
        point = []
        for coordinate in centroid:
            point.append(coordinate + random.gauss(0, 1))

        print(*point, sep=",", file=buf)

    buf.write("\n")
    buf.flush()
    buf.detach()
    file.seek(0)  # Return reading position to the start of the file
    return file


def make_dummy_input_file():
    file = tempfile.TemporaryFile()
    buf = io.TextIOWrapper(file)
    buf.write(DUMMY_INPUT)
    buf.flush()
    buf.detach()
    file.seek(0)  # Return reading position to the start of the file
    return file


def run_test(tests_dir: Path, config: dict[str, str], test_type: TestType):
    success = True

    input_path = tests_dir / Path(f"input_{config['idx']}.txt")
    output_path = tests_dir / Path(f"output_{config['idx']}.txt")
    with input_path.open() as input_file:
        result = execute(config, test_type, input_file)

    if result.returncode:
        success = False
        print(f"process returned with code {result.returncode}")

    if result.stderr:
        success = False
        print("process had non-empty stderr:")
        print_white_on_red(result.stderr)

    # Compare outputs
    if output_path.read_text().rstrip() != result.stdout.rstrip():
        success = False
        print("process had mismatching output")

    if success:
        print_green("success")
    else:
        print_red("failure")


def run_negative_test(
    config: dict[str, str],
    test_type: TestType,
    expect_msg: Optional[str] = None,
    input_file=None,
):
    expected_failure = True

    result = execute(config, test_type, input_file)

    if result.returncode != 1:
        expected_failure = False
        print(f"process returned with code {result.returncode}, expected 1")

    if expect_msg is not None and result.stdout.rstrip("\n") != expect_msg:
        expected_failure = False
        print(f'incorrect error msg - should print "{expect_msg}", printed:')
        print_white_on_blue(result.stdout)

    if result.stderr:
        print("process had non-empty stderr:")
        print_white_on_red(result.stderr)

    if expected_failure:
        print_green("success: expected failure")
    else:
        print_red("failure: unexpected result")


def execute(config, test_type, input_file=None) -> subprocess.CompletedProcess[str]:
    match test_type:
        case TestType.PYTHON:
            args = ["python3", "kmeans.py", config["k"]]
        case TestType.C:
            args = ["./kmeans", config["k"]]

    if config.get("max_iter"):
        args.append(config["max_iter"])

    result = subprocess.run(
        args,
        stdin=input_file,
        capture_output=True,
        text=True,
    )

    return result


def main():
    args = setup_argparser().parse_args()

    tests_dir = Path(args.tests_dir)

    readme_path = tests_dir / Path("test_readme.txt")
    configs = []
    with readme_path.open() as f:
        for line in f:
            match = LINE_CONFIG_REGEX.match(line)
            if match:
                config = match.groupdict()
                configs.append(config)

    test_types = []
    if args.test_python:
        test_types.append(TestType.PYTHON)

    if args.test_c:
        test_types.append(TestType.C)

    # Execute "positive tests" (expect success)
    for config in configs:
        print("Test", config["idx"])
        for test_type in test_types:
            print("testing", "python" if test_type is TestType.PYTHON else "C")
            run_test(tests_dir, config, test_type)

        print()

    random_input = generate_random_input()
    dummy_input = make_dummy_input_file()

    # Test for equality between Python and C on generated data
    if args.test_python and args.test_c:
        print("Testing for equality of Python and C outputs")
        for label, input_data in (("random", random_input), ("dummy", dummy_input)):
            python_run = execute({"k": "3"}, TestType.PYTHON, input_data)
            input_data.seek(0)  # Restore file pointer to the beginning of the file
            c_run = execute({"k": "3"}, TestType.C, input_data)
            input_data.seek(0)

            skip = False

            if python_run.returncode or python_run.stderr:
                print_red(f"execution of Python code on {label} data failed")
                if python_run.returncode:
                    print(
                        f"process had return code {python_run.returncode}, expected 0"
                    )

                if python_run.stderr:
                    print("process had non-empty stderr:")
                    print_white_on_red(python_run.stderr)
                skip = True

            if c_run.returncode or c_run.stderr:
                print_red(f"execution of C code on {label} data failed")
                if c_run.returncode:
                    print(f"process had return code {c_run.returncode}, expected 0")

                if c_run.stderr:
                    print("process had non-empty stderr:")
                    print_white_on_red(c_run.stderr)
                skip = True

            if skip:
                continue

            if python_run.stdout == c_run.stdout:
                print_green(
                    f"success: both Python and C had the same output on {label} data"
                )
            else:
                print_red(
                    f"failure: Python and C had different outputs on {label} data"
                )
        print()

    random_input.close()

    # Test for integer paramters written as floats
    print("Testing for integer parameters written as floats (K=3.0, iter=42.0000)")
    _config = {"k": "3.0", "max_iter": "42.0000"}
    for test_type, applicable in (
        (TestType.PYTHON, args.test_python),
        (TestType.C, args.test_c),
    ):
        if not applicable:
            continue

        success = True

        print("Testing", "python" if test_type == TestType.PYTHON else "C")
        result = execute(_config, TestType.PYTHON, dummy_input)
        dummy_input.seek(0)

        if result.returncode:
            success = False
            print(f"process returned with error code {result.returncode}, expected 0")

        if result.stderr:
            success = False
            print("process had non-empty stderr:")
            print_white_on_red(result.stderr)

        if success:
            print_green("success")
        else:
            print_red("failure")

        print()

    # Execute "negative tests" (expect failure)

    print("Testing invalid parameters for `iter`")
    for invalid_param in INVALID_PARAMETERS_VALUES:
        print(f"K=3, iter={invalid_param}")
        for test_type, applicable in (
            (TestType.PYTHON, args.test_python),
            (TestType.C, args.test_c),
        ):
            if not applicable:
                continue

            print("Testing", "python" if test_type == TestType.PYTHON else "C")
            run_negative_test(
                {"k": "3", "max_iter": invalid_param},
                test_type,
                "Incorrect maximum iteration!",
                dummy_input,
            )
            dummy_input.seek(0)
    print()

    print("Testing invalid parameters for K")
    for invalid_param in INVALID_PARAMETERS_VALUES:
        print(f"K={invalid_param}")
        for test_type, applicable in (
            (TestType.PYTHON, args.test_python),
            (TestType.C, args.test_c),
        ):
            if not applicable:
                continue

            print("Testing", "python" if test_type == TestType.PYTHON else "C")
            run_negative_test(
                {"k": invalid_param},
                test_type,
                "Incorrect number of clusters!",
                dummy_input,
            )
            dummy_input.seek(0)

    print("K=42, N=4")
    for test_type, applicable in (
        (TestType.PYTHON, args.test_python),
        (TestType.C, args.test_c),
    ):
        if not applicable:
            continue

        print("Testing", "python" if test_type == TestType.PYTHON else "C")
        run_negative_test(
            {"k": "42"},
            test_type,
            "Incorrect number of clusters!",
            dummy_input,
        )
        dummy_input.seek(0)

    dummy_input.close()


if __name__ == "__main__":
    main()
