

class CovigatorException(Exception):
    pass


class CovigatorDatabaseConnectionException(CovigatorException):
    pass


class CovigatorPipelineError(CovigatorException):
    pass


class CovigatorDashBoardInitialisationError(CovigatorException):
    pass


class CovigatorExcludedAssemblySequence(CovigatorException):
    pass


class CovigatorNotSupportedVariant(CovigatorException):
    pass


class CovigatorQueryException(CovigatorException):
    pass


class CovigatorExcludedSampleTooEarlyDateException(CovigatorException):
    pass


class CovigatorExcludedSampleTooManyMutations(CovigatorException):
    pass


class CovigatorErrorProcessingCoverageResults(CovigatorException):
    pass


class CovigatorExcludedSampleNarrowCoverage(CovigatorException):
    pass


class CovigatorExcludedSampleBadQualityReads(CovigatorException):
    pass


class CovigatorDashboardMissingPrecomputedData(CovigatorException):
    pass