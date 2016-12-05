import unittest
from bio_gtk_programme import LnFrame, seq_len

class TestSL(unittest.Testcase):
	def test_seq_len(self):
		self.assertTrue(seq_len >0)

if __name__ == '__main__':
	unittest.main()
