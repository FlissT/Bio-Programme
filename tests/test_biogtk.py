import unittest
#from bio_gtk_programme import LnFrame, seq_len

class TestStringMethods(unittest.TestCase):
    def test_upper(self):
        self.assertEqual('foo'.upper(), 'FOO')
    def test_upper(self):
        self.assertTrue('FOO'.isupper())

if __name__ == '__main__ ':
    unittest.main()

        
