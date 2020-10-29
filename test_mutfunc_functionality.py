import unittest
import pandas as pd
import mutfunc_functionality


class TestMutfunc(unittest.TestCase):

    def test_extract_snps(self):
        df = pd.read_excel("test/test_mutation_file.xlsx", header=4, keep_default_na=False)
        result_df = mutfunc_functionality.extract_snps(df)
        expected_df = pd.read_csv("test/test_extract_snps.csv")
        assert len(expected_df) == len(result_df)

if __name__ == "__main__":
    unittest.main()