

class CovigatorDatabaseConnectionException(Exception):
    pass


class CovigatorPipelineError(Exception):
    pass


class CovigatorDashBoardInitialisationError(Exception):
    pass


class CovigatorExcludedAssemblySequence(Exception):
    pass


class CovigatorNotSupportedVariant(Exception):
    pass


class CovigatorQueryException(Exception):
    pass


class CovigatorExcludedSampleTooEarlyDateException(Exception):
    pass


class CovigatorExcludedSampleTooManyMutations(Exception):
    pass


class CovigatorErrorProcessingCoverageResults(Exception):
    pass


class CovigatorExcludedSampleNarrowCoverage(Exception):
    pass


class CovigatorExcludedSampleBadQualityReads(Exception):
    pass


class CovigatorDashboardMissingPrecomputedData(Exception):
    pass