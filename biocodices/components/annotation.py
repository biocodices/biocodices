class Annotation:
    def __init__(self, rs):
        self.rs = rs

    def __repr__(self):
        return '<{} for {}>'.format(self.__class__.__name__, self.rs)
