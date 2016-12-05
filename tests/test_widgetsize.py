import unittest
from bio_gtk_programme import MyWindow


class MyWindowTestCase(unittest.TestCase):
    def setUp(self):
        self.widget = MyWindow("the widget")


    def test_default_widget_size(self):
        self.assertEqual(self.widget.size(),(300,300), "incorrect default size")

    def test_widget_resize(self):
        self.widget.resize(100, 150)
        self.assertEqual(self.widget.size(), (100, 150), "wrong size after resize")



if __name__ == '__main__':
    unittest.main()
