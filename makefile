test: install
	pytest test

install:
	python setup.py build
	pip install . -U

annotate:
	cython -a pyteomics/*.pyx

