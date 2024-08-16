import unittest

class Testing(unittest.TestCase):
    def test_string(self):
        a = 'some'
        b = 'some'
        self.assertEqual(a, b)

    def test_list(self):
        a = [1, 2, 3]
        b = [1, 2, 3]
        self.assertEqual(a, b)

    def test_float(self):
        a = 0.5
        b = 0.5
        self.assertEqual(a, b)


    if __name__ == '__main__':
        unittest.main()