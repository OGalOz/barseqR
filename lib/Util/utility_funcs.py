import logging

def stop(stop_str):
    if isinstance(stop_str, str):
        raise Exception(stop_str)
    elif isinstance(stop_str, int):
        raise Exception(str(stop_str))
    else:
        raise Exception("Stopped. Unknown")

def printl(print_str):
    logging.info(print_str)
