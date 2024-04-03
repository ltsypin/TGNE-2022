import os

file_dir = os.path.dirname(os.path.abspath(__file__))

# hardware utils

def get_cpu_cores():
    # If you're using Linux or macOS
    if os.name == 'posix':
        return os.cpu_count()

    # If you're using Windows
    elif os.name == 'nt':
        return multiprocessing.cpu_count()

    # If the operating system is not recognized
    else:
        return "Unable to determine the number of CPU cores."

