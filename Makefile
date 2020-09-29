PYTHON=python

.PHONY: test clean

test:
	${PYTHON} submit_tests.py
	# ${PYTHON} -m unittest discover -s test

clean:
	rm -rf test_*.{err,out}