PYTHON=python

.PHONY: test clean cancel

test:
	${PYTHON} submit_tests.py
	# ${PYTHON} -m unittest discover -s test

clean:
	rm -rf test_*.{err,out}

cancel:
	scancel -n run_pyscreener.batch
