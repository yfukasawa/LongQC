import sys, os
from subprocess import Popen, PIPE

from logging import getLogger
logger = getLogger(__name__)

class LqExec:

    def __init__(self, bin_path, ncpu=1):
        self.bin_path = bin_path
        self.proc = None

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
            _arg = [self.bin_path] + list(args)
            command = " ".join([str(i) for i in list(args)])
            logger.info("below command is executed: %s" % command)
            logger.info("%s is started." % self.bin_path)
        except OSError as e:
            logger.error(e)
            if out:
                fout.close()
            if err:
                ferr.close()
        else:
            if out:
                fout.close()
            if err:
                ferr.close()

    # still buggy
    def communicate(self, *args, inp=None):
        try:
            self.proc = Popen([self.bin_path] + list(args), stdout=PIPE, stderr=PIPE, stdin=PIPE)
            (stdout, stderr) = self.proc.communicate(input=inp)
            _arg = [self.bin_path] + list(args)
            command = " ".join([str(i) for i in list(args)])
            logger.info("below command is executed: %s" % command)
            logger.info("%s is started." % self.bin_path)
            print(stdout, stderr)
        except OSError as e:
            logger.error(e)
        else:
            pass

    def set_stdin(self, string):
        if not self.proc:
            return None
        self.proc.stdin.write(string)

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
        return str(self.proc.pid)

    def get_bin_path(self):
        return self.bin_path


# stand alone
if __name__ == "__main__":
    pass
    le = LqExec("/usr/bin/ls")
    le.exec("-l", "-a", out="./exec.txt")
