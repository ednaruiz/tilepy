import logging as lg

__all__ = ["settup_logging", "DEFAULT_LOGGER"]

DEFAULT_LOGGER = lg.getLogger()

def setup_logging(options=None,
                  output_filename=None,
                  *,
                  filename=None,
                  to_stdout="INFO",
                  to_file="INFO",
                  format="%(levelname)-8s: %(message)s",
                  datefmt=None,
                  logger_levels={},
                  dict_config=None,
                  capture_warnings=None,
                  skip_setup=None):
    """Configure the :mod:`logging` module.

    The default logging setup is given by the following equivalent `dict_config`
    (here in [yaml]_ format for better readability).

    ..
        If you change the code block below, please also change the corresponding block
        in :doc:`/intro/logging`.

    .. code-block :: yaml

        version: 1  # mandatory for logging config
        disable_existing_loggers: False  # keep module-based loggers already defined!
        formatters:
            custom:
                format: "%(levelname)-8s: %(message)s"   # options['format']
        handlers:
            to_stdout:
                class: logging.StreamHandler
                level: INFO         # options['to_stdout']
                formatter: custom
                stream: ext://sys.stdout
            to_file:
                class: logging.FileHandler
                level: INFO         # options['to_file']
                formatter: custom
                filename: output_filename.log   # options['filename']
                mode: a
        root:
            handlers: [to_stdout, to_file]
            level: DEBUG

    .. note ::
        We **remove** any previously configured logging handlers.
        This is to handle the case when this function is called multiple times,
        e.g., because you run multiple :class:`~tenpy.simulations.simulation.Simulation`
        classes sequentially (e.g., :func:`~tenpy.simulations.simulation.run_seq_simulations`).

    .. deprecated :: 0.9.0
        The arguments were previously collected in a dicitonary `options`.
        Now they should be given directly as keyword arguments.

    Parameters
    ----------
    **kwargs :
        Keyword arguments as described in the options below.
    output_filename : None | str
        The filename for where results are saved. The :cfg:option:`log.filename` for the
        log-file defaults to this, but replacing the extension with ``.log``.

    Options
    -------
    .. cfg:config :: log

        skip_setup: bool
            If True, don't change anything in the logging setup; just return.
            This is usefull for testing purposes, where `pytest` handles the logging setup.
            All other options are ignored in this case.
        to_stdout : None | ``"DEBUG" | "INFO" | "WARNING" | "ERROR" | "CRITICAL"``
            If not None, print log with (at least) the given level to stdout.
        to_file : None | ``"DEBUG" | "INFO" | "WARNING" | "ERROR" | "CRITICAL"``
            If not None, save log with (at least) the given level to a file.
            The filename is given by `filename`.
        filename : str
            Filename for the logfile.
            It defaults  to `output_filename` with the extension replaced to ".log".
            If ``None``, no log-file will be created, even with `to_file` set.
        logger_levels : dict(str, str)
            Set levels for certain loggers, e.g. ``{'tenpy.tools.params': 'WARNING'}`` to suppress
            the parameter readouts logs.
            The keys of this dictionary are logger names, which follow the module structure in
            tenpy.
            For example, setting the level for `tenpy.simulations` will change the level
            for all loggers in any of those submodules, including the one provided as
            ``Simluation.logger`` class attribute. Hence, all messages from Simulation class
            methods calling ``self.logger.info(...)`` will be affected by that.
        format : str
            Formatting string, `fmt` argument of :class:`logging.Formatter`.
            You can for example use ``"{loglevel:.4s} {asctime} {message}"`` to include the time
            stamp of each message into the log - this is usefull to get an idea where code hangs.
            Find
            `allowed keys <https://docs.python.org/3/library/logging.html#logrecord-attributes>`_
            here. The style of the formatter is chosen depending on whether the format string
            contains ``'%' '{' '$'``, respectively.
        datefmt : str
            Formatting string for the `asctime` key in the `format`, e.g. ``"%Y-%m-%d %H:%M:%S"``,
            see :meth:`logging.Formatter.formatTime`.
        dict_config : dict
            Alternatively, a full configuration dictionary for :func:`logging.config.dictConfig`.
            If used, all other options except `skip_setup` and `capture_warnings` are ignored.
        capture_warnings : bool
            Whether to call :func:`logging.captureWarnings` to include the warnings into the log.
    """
    import logging.config
    if options is not None:
        warnings.warn("Give logging parameters directly as keyword arguments!", FutureWarning, 2)
        locals().update(**options)
    if filename is None:
        if output_filename is not None:
            root, ext = os.path.splitext(output_filename)
            assert ext != '.log'
            filename = root + '.log'
    if capture_warnings is None:
        capture_warnings = dict_config is not None or to_stdout or to_file
    if skip_setup is None:
        skip_setup = skip_logging_setup
    if skip_setup:
        return
    if dict_config is None:
        handlers = {}
        if to_stdout:
            handlers['to_stdout'] = {
                'class': 'logging.StreamHandler',
                'level': to_stdout,
                'formatter': 'custom',
                'stream': 'ext://sys.stdout',
            }
        if to_file and filename is not None:
            handlers['to_file'] = {
                'class': 'logging.FileHandler',
                'level': to_file,
                'formatter': 'custom',
                'filename': filename,
                'mode': 'a',
            }
            if not to_stdout:
                cwd = os.getcwd()
                print(f"now logging to {cwd!s}/{filename!s}")
        dict_config = {
            'version': 1,  # mandatory
            'disable_existing_loggers': False,
            'formatters': {
                'custom': {
                    'format': format,
                    'datefmt': datefmt
                }
            },
            'handlers': handlers,
            'root': {
                'handlers': list(handlers.keys()),
                'level': 'DEBUG'
            },
            'loggers': {},
        }
        if '%' not in format:
            if '{' in format:
                assert '$' not in format
                style = '{'
            else:
                assert '$' in format
                style = '$'
            dict_config['formatters']['custom']['style'] = style
        for name, level in logger_levels.items():
            if name == 'root':
                dict_config['root']['level'] = level
            else:
                dict_config['loggers'].setdefault(name, {})['level'] = level
    else:
        dict_config.setdefault('disable_existing_loggers', False)
    # note: dictConfig cleans up previously existing handlers etc
    logging.config.dictConfig(dict_config)
    if capture_warnings:
        logging.captureWarnings(True)

