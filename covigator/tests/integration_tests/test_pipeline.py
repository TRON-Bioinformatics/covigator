import unittest
from covigator.processor.pipeline import Pipeline, CovigatorPipelineError
from covigator.processor.gisaid_pipeline import GisaidPipeline

class PipelineTest(unittest.TestCase):

    def test_pipeline_run(self):
        fastq1 = "SRR11140748_R1.fastq.gz"
        fastq2 = "SRR11140748_R2.fastq.gz"
        vcf_file = "expected_snpeff.vcf"
        p = Pipeline()
    
        # TODO: Implement Unit Test here
        # Maybe do diff as assertion?    
        #self.assertEqual(p.run(fastq1=fastq1, fastq2=fastq2), vcf_file)
        #with self.assertRaises(CovigatorPipelineError):
            #vcf_file = p.run(fastq1=fastq1, fastq2=fastq2)
        self.assertEqual(open(p.run(fastq1=fastq1, fastq2=fastq2)).read(), open(vcf_file).read())

class GisaidPipelineTest(unittest.TestCase):
    
    def test_pipeline_run(self):
        run_accession = "EPI_ISL_417140"
        vcf_file = "gisaid.vcf"

        p = GisaidPipeline()

        self.assertEqual(open(p.run(run_accession=run_accession)).read(), open(vcf_file).read())

if __name__ == '__main__':
    unittest.main()
