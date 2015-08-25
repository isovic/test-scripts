#! /usr/bin/python

def execute_command(command):
	sys.stderr.write('[%s wrapper] %s\n' % (MAPPER_NAME, command));
	subprocess.call(command, shell=True);
	sys.stderr.write('\n');

def execute_command_get_stdout(command):
	sys.stderr.write('[%s wrapper] %s\n' % (MAPPER_NAME, command));
	p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
	[out, err] = p.communicate()
	sys.stderr.write('\n');

	return [out, err];

def main():
	pass;

if __name__ == "__main__":
	main()
