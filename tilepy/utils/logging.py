"""Legacy module for the setup of logging for argparse-based scripts."""
import logging as lg

__all__ = ["setup_logging", "DEFAULT_LOGGER"]

DEFAULT_LOGGER = lg.getLogger()


def setup_logging(logger, log_path, console_log_level, overwrite=False):
    """Set up a logger.

    Parameters
    ----------
    logger: logging.Logger
        Default logger.
    log_path: pathlib.Path
        Output path for the log file.
    console_log_level: str
        Level for the logging to console.
        File records always everything.

    Returns
    -------
    logger: logging.Logger
        Initialized logger.
    """
    if (log_path).exists() and not overwrite:
        raise FileExistsError("Found previous log file. Move it or use --overwrite.")

    logger.setLevel(lg.DEBUG)
    numeric_level = getattr(lg, console_log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: %s" % console_log_level)
    # Create handlers
    console_handler = lg.StreamHandler()
    writing_mode = "w" if overwrite else "a"
    file_handler = lg.FileHandler(log_path, mode=writing_mode)
    console_handler.setLevel(numeric_level)
    file_handler.setLevel(lg.DEBUG)
    # Create formatters and add it to handlers
    c_format = lg.Formatter("%(module)s - %(funcName)s - %(levelname)s - %(message)s")
    f_format = lg.Formatter("%(module)s - %(funcName)s - %(levelname)s - %(message)s")
    console_handler.setFormatter(c_format)
    file_handler.setFormatter(f_format)
    # Add handlers to the logger
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

    return logger
