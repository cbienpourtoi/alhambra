__author__ = 'loic'

from sklearn.linear_model import LogisticRegression
from sklearn import preprocessing


class machine(object):
    def __init__(self):
        self.scaler = None
        self.model = None

    def train(self, X, Y):
        """
        :param X: training sample X
        :param Y: training sample classes
        :return:
        """
        raise NotImplementedError()


    def use(self, X):
        """
        :param X: unknown samples
        :return: predicted classes
        """
        Xnorm = self.scaler.transform(X)
        y_prediction = self.model.predict(Xnorm)
        return y_prediction

    def _train_scaler(self, X):
        self.scaler = preprocessing.StandardScaler().fit(X)


class logistic(machine):
    def train(self, X, Y):
        self._train_scaler(X)
        Xnorm = self.scaler.transform(X)
        self.model = LogisticRegression(C=1.)
        self.model.fit(Xnorm, Y)


class lda(machine):
    def train(self, X, Y):
        self._train_scaler(X)
        Xnorm = self.scaler.transform(X)
        self.model = LDA(C=1.)
        self.model.fit(Xnorm, Y)
