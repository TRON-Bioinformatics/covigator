

class CovigatorException(Exception):
    pass


class CovigatorDatabaseConnectionException(CovigatorException):
    pass


class CovigatorPipelineError(CovigatorException):
    pass


class CovigatorDashBoardInitialisationError(CovigatorException):
    pass


class CovigatorNotSupportedVariant(CovigatorException):
    pass


class CovigatorQueryException(CovigatorException):
    pass


class CovigatorErrorProcessingCoverageResults(CovigatorException):
    pass


class CovigatorDashboardMissingPrecomputedData(CovigatorException):
    pass


class CovigatorExcludedSampleException(CovigatorException):
    pass


class CovigatorExcludedAssemblySequence(CovigatorExcludedSampleException):
    pass


class CovigatorExcludedSampleTooEarlyDateException(CovigatorExcludedSampleException):
    pass


class CovigatorExcludedSampleTooManyMutations(CovigatorExcludedSampleException):
    pass


class CovigatorExcludedSampleNarrowCoverage(CovigatorExcludedSampleException):
    pass


class CovigatorExcludedSampleBadQualityReads(CovigatorExcludedSampleException):
    pass
