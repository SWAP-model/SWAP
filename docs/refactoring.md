# SWAP refactoring

When trying to crate SWAP models that smoothly and efficiently run in parallel, it turned out that the SWAP code design follows too many legacy patterns that inhibit full parallelization. Lack of accessible documentation also caused 