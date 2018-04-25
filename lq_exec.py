import sys, os
from subprocess import Popen, PIPE

class LqExec:

    def __init__(self, bin_path, logger=None):
        self.bin_path = bin_path
        self.proc = None
        if logger:
            self.logger = logger
        else:
            self.logger = logging.getLogger(__name__)
            logger.setLevel(logging.INFO)
            sh = logging.StreamHandler()
            formatter = logging.Formatter('%(module)s:%(asctime)s:%(lineno)d:%(levelname)s:%(message)s')
            sh.setFormatter(formatter)
            logger.addHandler(sh)

    def exec(self, *args, out=None, err=None):
        if not out:
            fout = PIPE
        else:
            fout = open(out, "w")
        if not err:
            ferr = PIPE
        else:
            ferr = open(err, "w")

        try:
            self.proc = Popen([self.bin_path] + list(args), stdout=fout, stderr=ferr)
            if self.logger:
                _arg = [self.bin_path] + list(args)
                command = " ".join([str(i) for i in list(args)])
                self.logger.debug("below command is executed: %s" % command)
                self.logger.info("%s is started." % self.bin_path)
        except:
            if out:
                fout.close()
            if err:
                ferr.close()
        else:
            if out:
                fout.close()
            if err:
                ferr.close()

    def get_stdout(self):
        if not self.proc:
            return None
        return self.proc.stdout

    def get_stderr(self):
        if not self.proc:
            return None
        return self.proc.sterr

    def get_poll(self):
        return self.proc.poll()

    def get_pid(self):
        return self.proc.pid

    def get_bin_path(self):
        return self.bin_path


# stand alone
if __name__ == "__main__":

    #print(subprocess.check_output('ls'))

    le = LqExec("/usr/bin/ls")
    le.exec("-ll", "-a", out="/home/fukasay/test_exec.txt")



