import argparse
from datetime import datetime, timedelta
from collections import defaultdict

from tabulate import tabulate

FMT = "%Y_%m_%d-%H:%M:%S"


def parse_times(lines):
    time_per_day = defaultdict(timedelta)

    prev_start = None

    for line in lines:
        time_str, what = line.split()
        time = datetime.strptime(time_str, FMT)

        date = time.date()

        if prev_start is None:
            assert what == "start"
            prev_start = time

        else:
            assert what == "end"
            time_per_day[date] += time - prev_start
            prev_start = None

    return time_per_day


def report(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    time_per_day = parse_times(lines)

    time_per_week = defaultdict(timedelta)

    for date, time in time_per_day.items():
        time_per_week[date.isocalendar().week] += time


    def hours(time):
        seconds = time.seconds
        h, remainder = divmod(time.seconds, 3600)
        m, s= divmod(remainder, 60)
        return f"{time.days*24 + h}:{m:02d}:{s:02d}"

    time_str_per_day = [(day, hours(time)) for day, time in time_per_day.items()]
    time_str_per_week = [(week, hours(time)) for week, time in time_per_week.items()]

    print(tabulate(time_str_per_day, headers=["day", "hours"]))
    print()
    print(tabulate(time_str_per_week, headers=["week", "hours"]))


def cli():
    parser = argparse.ArgumentParser(
        prog="Arbeitszeit",
        description="",
        epilog="TT 2023",
    )

    parser.add_argument(
        "what", type=str, nargs="?", default="report"
    )  # positional argument
    parser.add_argument("-f", "--file", type=str, default="Arbeitszeit.dat")

    args = parser.parse_args()

    if args.what in ["start", "end"]:
        msg = f"{datetime.now().strftime(FMT)} {args.what}"
        with open(args.file, "a") as f:
            f.write(msg + "\n")
        print(f"{msg} --> {args.file}")

    report(args.file)


if __name__ == "__main__":
    cli()
