from threading import Thread
from read_evaluator import evaluate_pair
__author__ = 'jgrundst'


class EvaluatorWorker(Thread):

    def __init__(self, paired=True, queue=None, FW=None):
        """
        Threaded worker class to perform read evaluations
        assumes paired end.
        :param paired:
        :param queue:
        :param FW: FastqWriter from the BamParser class
        :return:
        """
        Thread.__init__(self)
        self.queue = queue
        self.FW = FW

    def run(self):
        while True:
            read1, read2 = self.queue.get()
            if evaluate_pair(read1, read2):
                self.FW.print_reads(read1=read1, read2=read2)
