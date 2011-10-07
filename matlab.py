import subprocess


startupcmds = ['cd /Applications/tomlab/',
               'startup',
               'cd ~/mitnotes/research/provenance/src/sim',
               'magic(10)']

args = ['/Applications/MATLAB_R2010b.app/bin/matlab', '-nojvm']
p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
#print p.communicate()

print "===="
print p.poll()

for cmd in startupcmds:
    print "command", p.stdin.write('%s\n' % cmd)


print "reading from pipe"
p.communicate('exit\n')
p.terminate()
