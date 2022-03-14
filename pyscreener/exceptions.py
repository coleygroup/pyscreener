class BadPDBFileError(Exception):
    pass


class MissingExecutableError(Exception):
    pass


class MissingEnvironmentVariableError(Exception):
    pass


class MissingFileError(Exception):
    pass


class MisconfiguredDirectoryError(Exception):
    pass


class NotSimulatedError(Exception):
    pass


class InvalidResultError(Exception):
    pass


class UnsupportedSoftwareError(Exception):
    pass


class ReceptorPreparationError(Exception):
    pass
