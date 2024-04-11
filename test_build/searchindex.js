Search.setIndex({"docnames": ["api/executor", "api/fields", "api/index", "contributing", "index", "installation"], "filenames": ["api/executor.rst", "api/fields.rst", "api/index.rst", "contributing.rst", "index.rst", "installation.rst"], "titles": ["Executor", "Fields", "API", "Contributing", "Welcome to NeoFOAM!", "Installation"], "terms": {"neofoam": [0, 1, 3, 5], "us": [0, 1, 3, 5], "mpi": 0, "x": 0, "approach": 0, "parallel": 0, "where": 0, "i": [0, 1, 3, 4], "execut": [0, 3, 5], "space": 0, "The": [0, 1, 3, 4, 5], "class": [0, 1], "kokko": [0, 1, 5], "provid": [0, 5], "an": [0, 1], "interfac": [0, 2], "memori": [0, 1], "manag": 0, "specif": 0, "were": 0, "oper": [0, 1], "cpuexecutor": 0, "run": 0, "cpu": [0, 1], "ompexecutor": 0, "openmp": [0, 1, 5], "gpuexecutor": [0, 1], "gpu": [0, 1, 4], "One": 0, "goal": [0, 4], "abil": 0, "quickli": 0, "switch": [0, 5], "between": [0, 5], "model": [0, 1], "runtim": 0, "gpuexec": [0, 1], "cpuexec": 0, "field": [0, 2, 4], "scalar": [0, 1], "gpufield": 0, "10": 0, "cpufield": 0, "std": [0, 1], "variant": 0, "allow": [0, 1], "differ": [0, 5], "strategi": 0, "alloc": [0, 1], "we": [0, 1, 4, 5], "visit": 0, "exec": [0, 1], "const": [0, 1], "auto": [0, 1], "functor": [0, 1], "ar": [0, 1, 3, 4, 5], "struct": 0, "void": [0, 1], "cout": 0, "endl": 0, "pattern": 0, "abov": 0, "would": 0, "print": 0, "messag": 0, "depend": [0, 5], "type": 0, "To": [0, 5], "extend": 0, "librari": [0, 4], "addit": [0, 1], "featur": 0, "should": [0, 1, 5], "implement": [0, 1, 4], "can": [0, 3, 5], "check": 0, "two": [0, 1], "same": [0, 1], "e": [0, 4], "equal": 0, "central": 1, "element": 1, "platform": [1, 4], "portabl": [1, 4], "cfd": 1, "framework": 1, "perform": 1, "basic": [1, 4], "algebra": 1, "binari": [1, 4], "like": [1, 4], "subtract": 1, "multipl": 1, "In": 1, "follow": [1, 3, 5], "explain": 1, "detail": [1, 5], "exampl": [1, 4], "block": 1, "code": [1, 4], "below": [1, 5], "show": 1, "nodiscard": 1, "t": [1, 4], "rh": 1, "result": 1, "exec_": 1, "size_": 1, "thi": [1, 4, 5], "add": [1, 4], "return": 1, "besid": 1, "creat": 1, "temporari": 1, "mainli": 1, "call": 1, "free": 1, "stand": 1, "function": 1, "which": [1, 4], "fieldoper": 1, "hpp": 1, "turn": 1, "dispatch": 1, "addop": 1, "hold": 1, "actual": 1, "kernel": 1, "case": 1, "parallel_for": 1, "see": 1, "document": [1, 4, 5], "more": 1, "executor": [1, 2, 4], "typenam": 1, "a_f": 1, "b_f": 1, "b": [1, 3], "rangepolici": 1, "0": 1, "size": 1, "kokkos_class_lambda": 1, "int": 1, "snippet": 1, "also": 1, "highlight": 1, "anoth": 1, "import": 1, "aspect": 1, "here": [1, 3], "defin": 1, "program": 1, "deallocationg": 1, "A": 1, "full": 1, "could": [1, 4], "gpua": 1, "n": 1, "fill": 1, "1": 1, "gpub": 1, "2": 1, "gpuc": 1, "namespac": 1, "templat": 1, "includ": 1, "contain": 1, "data": 1, "some": [1, 3], "public": 1, "inlin": 1, "size_t": 1, "given": 1, "paramet": 1, "associ": 1, "matrix": 1, "destroi": 1, "object": 1, "func": 1, "appli": 1, "f": 1, "transform": 1, "ideal": 1, "kokkos_lamba": 1, "map": 1, "over": 1, "copytohost": 1, "copi": 1, "back": 1, "host": 1, "from": [1, 4, 5], "anywher": 1, "pars": 1, "exit": 1, "sourc": [1, 4], "must": 1, "kokkos_funct": 1, "index": [1, 3, 4], "cell": 1, "valu": [1, 5], "assign": 1, "set": [1, 4, 5], "arithmet": 1, "second": [1, 3], "multipli": 1, "everi": 1, "setsiz": 1, "new": 1, "direct": 1, "access": 1, "underli": 1, "pointer": 1, "first": [1, 3], "get": [1, 3], "refer": 1, "span": 1, "privat": 1, "member": 1, "data_": 1, "etc": 1, "overview": 2, "design": 2, "highli": [3, 4], "welcom": 3, "inform": 3, "you": [3, 4, 5], "start": 3, "found": 3, "main": 3, "doxygen": 3, "onlin": 3, "howev": [3, 4, 5], "want": [3, 4], "local": 3, "do": 3, "so": [3, 4], "step": [3, 5], "make": 3, "sure": 3, "sphinx": 3, "instal": [3, 4], "your": 3, "system": 3, "command": [3, 5], "cmake": [3, 4], "dneofoam_build_doc": 3, "ON": [3, 5], "configur": [3, 5], "target": 3, "built": 3, "doc": 3, "directori": [3, 5], "view": 3, "open": [3, 5], "html": 3, "file": [3, 4], "web": 3, "browser": 3, "firefox": 3, "altern": 3, "just": 3, "ad": 3, "process": 3, "project": [4, 5], "ha": 4, "itself": 4, "bring": 4, "modern": 4, "softwar": 4, "develop": 4, "method": 4, "core": 4, "By": [4, 5], "reimplement": 4, "libfinitevolum": 4, "libopenfoam": 4, "deliv": 4, "compliant": 4, "c": 4, "20": 4, "extens": 4, "unit": 4, "test": [4, 5], "readi": 4, "via": 4, "plugin": 4, "aim": 4, "high": 4, "level": 4, "interoper": 4, "reason": 4, "might": 4, "deviat": 4, "api": 4, "commun": 4, "driven": 4, "contribut": 4, "everyon": 4, "preset": 4, "don": 4, "expect": 4, "abi": 4, "mean": 4, "won": 4, "produc": 4, "serv": 4, "replac": 4, "exist": [4, 5], "instead": 4, "possibl": 4, "compil": [4, 5], "pimplefoam": 4, "other": 4, "against": 4, "current": 4, "support": [4, 5], "veri": 4, "limit": 4, "simplest": 4, "wai": [4, 5], "adopt": 4, "procedur": [4, 5], "showcas": 4, "our": 4, "folder": 4, "For": 4, "now": 4, "onli": 4, "have": [4, 5], "minim": 4, "still": 4, "work": [4, 5], "But": 4, "gradual": 4, "move": 4, "along": 4, "most": [4, 5], "initi": 4, "port": 4, "solver": 4, "sinc": 4, "pre": 4, "postprocess": 4, "tool": 4, "basi": 4, "trivial": 4, "modul": 4, "search": 4, "page": 4, "clone": 5, "repositori": 5, "git": 5, "http": 5, "github": 5, "com": 5, "exasim": 5, "navig": 5, "cd": 5, "thu": 5, "standard": 5, "recommend": 5, "one": 5, "mkdir": 5, "desiredbuildflag": 5, "chain": 5, "d": 5, "relev": 5, "flag": 5, "descript": 5, "default": 5, "cmake_build_typ": 5, "debug": 5, "releas": 5, "mode": 5, "neofoam_build_app": 5, "applic": 5, "neofoam_build_benchmark": 5, "benchmark": 5, "off": 5, "neofoam_build_doc": 5, "neofoam_build_test": 5, "kokkos_enable_seri": 5, "enabl": 5, "serial": 5, "backend": 5, "kokkos_enable_openmp": 5, "kokkos_enable_rocm": 5, "rocm": 5, "kokkos_enable_sycl": 5, "sycl": 5, "kokkos_enable_cuda": 5, "cuda": 5, "gui": 5, "easili": 5, "addition": 5, "sever": 5, "commonli": 5, "requir": 5, "combin": 5, "list": 5, "ninja": 5, "common": 5, "devic": 5, "It": 5, "note": 5, "chang": 5, "chosen": 5, "them": 5}, "objects": {"": [[1, 0, 1, "_CPPv47NeoFOAM", "NeoFOAM"], [1, 1, 1, "_CPPv4I0EN7NeoFOAM5FieldE", "NeoFOAM::Field"], [1, 2, 1, "_CPPv4N7NeoFOAM5Field5FieldERK8executor6size_t", "NeoFOAM::Field::Field"], [1, 3, 1, "_CPPv4N7NeoFOAM5Field5FieldERK8executor6size_t", "NeoFOAM::Field::Field::exec"], [1, 3, 1, "_CPPv4N7NeoFOAM5Field5FieldERK8executor6size_t", "NeoFOAM::Field::Field::size"], [1, 4, 1, "_CPPv4I0EN7NeoFOAM5FieldE", "NeoFOAM::Field::T"], [1, 2, 1, "_CPPv4I0EN7NeoFOAM5Field5applyEv4func", "NeoFOAM::Field::apply"], [1, 3, 1, "_CPPv4I0EN7NeoFOAM5Field5applyEv4func", "NeoFOAM::Field::apply::f"], [1, 4, 1, "_CPPv4I0EN7NeoFOAM5Field5applyEv4func", "NeoFOAM::Field::apply::func"], [1, 2, 1, "_CPPv4N7NeoFOAM5Field10copyToHostER5FieldI1TE", "NeoFOAM::Field::copyToHost"], [1, 2, 1, "_CPPv4N7NeoFOAM5Field10copyToHostEv", "NeoFOAM::Field::copyToHost"], [1, 3, 1, "_CPPv4N7NeoFOAM5Field10copyToHostER5FieldI1TE", "NeoFOAM::Field::copyToHost::result"], [1, 2, 1, "_CPPv4N7NeoFOAM5Field4dataEv", "NeoFOAM::Field::data"], [1, 2, 1, "_CPPv4NK7NeoFOAM5Field4dataEv", "NeoFOAM::Field::data"], [1, 5, 1, "_CPPv4N7NeoFOAM5Field5data_E", "NeoFOAM::Field::data_"], [1, 2, 1, "_CPPv4N7NeoFOAM5Field4execEv", "NeoFOAM::Field::exec"], [1, 5, 1, "_CPPv4N7NeoFOAM5Field5exec_E", "NeoFOAM::Field::exec_"], [1, 2, 1, "_CPPv4N7NeoFOAM5Field5fieldEv", "NeoFOAM::Field::field"], [1, 2, 1, "_CPPv4NK7NeoFOAM5Field5fieldEv", "NeoFOAM::Field::field"], [1, 2, 1, "_CPPv4N7NeoFOAM5FieldmlEK6scalar", "NeoFOAM::Field::operator*"], [1, 2, 1, "_CPPv4N7NeoFOAM5FieldmlERK5FieldI6scalarE", "NeoFOAM::Field::operator*"], [1, 3, 1, "_CPPv4N7NeoFOAM5FieldmlEK6scalar", "NeoFOAM::Field::operator*::rhs"], [1, 3, 1, "_CPPv4N7NeoFOAM5FieldmlERK5FieldI6scalarE", "NeoFOAM::Field::operator*::rhs"], [1, 2, 1, "_CPPv4N7NeoFOAM5FieldplERK5FieldI1TE", "NeoFOAM::Field::operator+"], [1, 3, 1, "_CPPv4N7NeoFOAM5FieldplERK5FieldI1TE", "NeoFOAM::Field::operator+::rhs"], [1, 2, 1, "_CPPv4N7NeoFOAM5FieldmiERK5FieldI1TE", "NeoFOAM::Field::operator-"], [1, 3, 1, "_CPPv4N7NeoFOAM5FieldmiERK5FieldI1TE", "NeoFOAM::Field::operator-::rhs"], [1, 2, 1, "_CPPv4N7NeoFOAM5FieldaSERK1T", "NeoFOAM::Field::operator="], [1, 2, 1, "_CPPv4N7NeoFOAM5FieldaSERK5FieldI1TE", "NeoFOAM::Field::operator="], [1, 3, 1, "_CPPv4N7NeoFOAM5FieldaSERK1T", "NeoFOAM::Field::operator=::rhs"], [1, 3, 1, "_CPPv4N7NeoFOAM5FieldaSERK5FieldI1TE", "NeoFOAM::Field::operator=::rhs"], [1, 2, 1, "_CPPv4N7NeoFOAM5Field7setSizeEK6size_t", "NeoFOAM::Field::setSize"], [1, 3, 1, "_CPPv4N7NeoFOAM5Field7setSizeEK6size_t", "NeoFOAM::Field::setSize::size"], [1, 2, 1, "_CPPv4NK7NeoFOAM5Field4sizeEv", "NeoFOAM::Field::size"], [1, 5, 1, "_CPPv4N7NeoFOAM5Field5size_E", "NeoFOAM::Field::size_"], [1, 2, 1, "_CPPv4N7NeoFOAM5FieldD0Ev", "NeoFOAM::Field::~Field"]]}, "objtypes": {"0": "cpp:type", "1": "cpp:class", "2": "cpp:function", "3": "cpp:functionParam", "4": "cpp:templateParam", "5": "cpp:member"}, "objnames": {"0": ["cpp", "type", "C++ type"], "1": ["cpp", "class", "C++ class"], "2": ["cpp", "function", "C++ function"], "3": ["cpp", "functionParam", "C++ function parameter"], "4": ["cpp", "templateParam", "C++ template parameter"], "5": ["cpp", "member", "C++ member"]}, "titleterms": {"executor": 0, "overview": [0, 1], "design": 0, "field": 1, "interfac": 1, "api": 2, "contribut": 3, "build": [3, 4, 5], "document": 3, "welcom": 4, "neofoam": 4, "tabl": 4, "content": 4, "compat": 4, "openfoam": 4, "applic": 4, "indic": 4, "instal": 5, "cmake": 5, "preset": 5}, "envversion": {"sphinx.domains.c": 3, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 9, "sphinx.domains.index": 1, "sphinx.domains.javascript": 3, "sphinx.domains.math": 2, "sphinx.domains.python": 4, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.intersphinx": 1, "sphinx.ext.todo": 2, "sphinx.ext.viewcode": 1, "sphinx": 60}, "alltitles": {"Executor": [[0, "executor"]], "Overview": [[0, "overview"], [1, "overview"]], "Design": [[0, "design"]], "Fields": [[1, "fields"]], "Interface": [[1, "interface"]], "API": [[2, "api"]], "Contributing": [[3, "contributing"]], "Building the Documentation": [[3, "building-the-documentation"]], "Welcome to NeoFOAM!": [[4, "welcome-to-neofoam"]], "Table of Contents": [[4, "table-of-contents"]], "Compatibility with OpenFOAM": [[4, "compatibility-with-openfoam"]], "Building OpenFOAM Applications with NeoFOAM": [[4, "building-openfoam-applications-with-neofoam"]], "Indices and tables": [[4, "indices-and-tables"]], "Installation": [[5, "installation"]], "Building with Cmake Presets": [[5, "building-with-cmake-presets"]]}, "indexentries": {"neofoam (c++ type)": [[1, "_CPPv47NeoFOAM"]], "neofoam::field (c++ class)": [[1, "_CPPv4I0EN7NeoFOAM5FieldE"]], "neofoam::field::field (c++ function)": [[1, "_CPPv4N7NeoFOAM5Field5FieldERK8executor6size_t"], [1, "_CPPv4N7NeoFOAM5Field5fieldEv"], [1, "_CPPv4NK7NeoFOAM5Field5fieldEv"]], "neofoam::field::apply (c++ function)": [[1, "_CPPv4I0EN7NeoFOAM5Field5applyEv4func"]], "neofoam::field::copytohost (c++ function)": [[1, "_CPPv4N7NeoFOAM5Field10copyToHostER5FieldI1TE"], [1, "_CPPv4N7NeoFOAM5Field10copyToHostEv"]], "neofoam::field::data (c++ function)": [[1, "_CPPv4N7NeoFOAM5Field4dataEv"], [1, "_CPPv4NK7NeoFOAM5Field4dataEv"]], "neofoam::field::data_ (c++ member)": [[1, "_CPPv4N7NeoFOAM5Field5data_E"]], "neofoam::field::exec (c++ function)": [[1, "_CPPv4N7NeoFOAM5Field4execEv"]], "neofoam::field::exec_ (c++ member)": [[1, "_CPPv4N7NeoFOAM5Field5exec_E"]], "neofoam::field::operator* (c++ function)": [[1, "_CPPv4N7NeoFOAM5FieldmlEK6scalar"], [1, "_CPPv4N7NeoFOAM5FieldmlERK5FieldI6scalarE"]], "neofoam::field::operator+ (c++ function)": [[1, "_CPPv4N7NeoFOAM5FieldplERK5FieldI1TE"]], "neofoam::field::operator- (c++ function)": [[1, "_CPPv4N7NeoFOAM5FieldmiERK5FieldI1TE"]], "neofoam::field::operator= (c++ function)": [[1, "_CPPv4N7NeoFOAM5FieldaSERK1T"], [1, "_CPPv4N7NeoFOAM5FieldaSERK5FieldI1TE"]], "neofoam::field::setsize (c++ function)": [[1, "_CPPv4N7NeoFOAM5Field7setSizeEK6size_t"]], "neofoam::field::size (c++ function)": [[1, "_CPPv4NK7NeoFOAM5Field4sizeEv"]], "neofoam::field::size_ (c++ member)": [[1, "_CPPv4N7NeoFOAM5Field5size_E"]], "neofoam::field::~field (c++ function)": [[1, "_CPPv4N7NeoFOAM5FieldD0Ev"]]}})