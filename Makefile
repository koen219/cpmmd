.PHONY : lib/pycpm python clean all scripts hoomd
all: scripts

PY=python3.12

PYTHON_VERSION = $(shell ${PY} -c 'import sys; print("python{}.{}".format(*sys.version_info[0:2]))')
VENV_PKG = .venv/lib/$(PYTHON_VERSION)/site-packages
VENV_HOOMD = $(VENV_PKG)/hoomd/__init__.py

$(VENV_HOOMD):
	. .venv/bin/activate
	cd lib/hoomd && make

hoomd: $(VENV_HOOMD)

scripts: lib/pycpm hoomd 
	.venv/bin/python -m pip install -e .

lib/pycpm: python .venv/bin/maturin 
	. .venv/bin/activate && cd lib/pycpm2 && maturin develop -r

.venv/bin/maturin: .venv/bin/activate
	. .venv/bin/activate
	.venv/bin/python -m pip install maturin
	.venv/bin/python -m pip install -e lib/ECMGen
	touch .venv/bin/maturin

.venv/bin/activate: requirements.txt
	test -d .venv || ${PY} -m venv ./.venv
	.venv/bin/python -m pip install pip --upgrade
	.venv/bin/python -m pip install -Ur requirements.txt
	touch .venv/bin/activate

clean:
	rm -rf .venv
	cd src/cpmmd && cargo clean
	cd lib/cpm2 && cargo clean
	cd lib/pycpm2 && cargo clean
