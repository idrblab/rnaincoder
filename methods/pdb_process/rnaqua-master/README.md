# A tool supporting quality assessment of RNA 3D structures

RNA QUality Assessment tool (RNAQUA) is RESTful web service client developed in Java using [Jersey][jersey] providing a set of web services designed initially for [RNAssess][rnassess] to support quality assessment of RNA 3D structures. At the input, RNA 3D structures must be introduced in PDB format. First, a comprehensive validation of input RNA 3D structure(s) is performed and all identified inconsistencies are reported to a user. An output is returned in XML format using [JAXB][jaxb]. Features of the tool include:

1. A single 3D structure-based analysis.
  - A comprehensive validation of the input RNA 3D structure(s) to report all PDB format inconsistencies (*PDB-VALIDATION*). 
  - A computation of an overall clash score (*CLASH-SCORE*), i.e., number of bad overlaps per 1000 atoms, provided by [MolProbity][molprobity].
  - A sequence-based analysis allowing to investigate a sequence of RNA 3D structure(s) (*SEQUENCE*) or a set of incontinuous RNA 3D substructures specified by the user (*FRAGMENT*).
  - An extraction (*ORIGINAL-3D*) and unification (*RENUMERATED-3D*) of RNA 3D structure(s) useful during analysis of input RNA 3D models or a set of incontinuous 3D substructures differing in sequence, distribution of chains, and numbering of residues.

2. Analysis of 3D model(s) within the reference structure context.
  - A computation of the following measures: [Root-mean-square deviation][rmsd] (*ROOT-MEAN-SQUARE-DEVIATION*), [Interaction network fidelity][inf] (*ALL-INTERACTION-NETWORK-FIDELITY-SCORES-AT-ONCE*,
*INTERACTION-NETWORK-FIDELITY-WATSON-CRICK*, *INTERACTION-NETWORK-FIDELITY-NON-WATSON-CRICK*, *INTERACTION-NETWORK-FIDELITY-STACKING*, *INTERACTION-NETWORK-FIDELITY-ALL*), [Deformation index][inf] (*DEFORMATION-INDEX*), [P-value][pvalue] (*P-VALUE*), all scores at once (*ALL-SCORES-AT-ONCE*).
  - A sequence-based comparison performed between input RNA 3D model(s) (*SEQUENCE*) or a set of incontinuous 3D substructures specified by the user (*FRAGMENT*) and the corresponding 3D structure/substructures of the reference to investigate disturbances in the compared sequences.
  - An extraction (*ORIGINAL-3D*) and unification (*RENUMERATED-3D*) of RNA 3D model(s) or a set of incontinuous 3D substructures specified by the user which are additionally superimposed over the corresponding 3D structure/substructures of the reference. At the output, ZIP archive including the coordinates of the reference structure as well as all considered RNA 3D model(s) is returned.   
  - A generation of [Deformation profile][inf] (*DEFORMATION-PROFILE*) provided by [DPS_1.0.1][inf] which was made available by Jose Almeida Cruz and Eric Westhof.
  
3. To ensure robustness of quality assessment process a user can specify the appropriate alignment (*-a,--alignment*) between the reference 3D structure and all analyzed RNA 3D model(s) which often differ slightly in sequence, distribution of chains or numbering of residues. An example of alignment prepared between the reference structure (*solution.pdb*) and a single RNA 3D model (*model.pdb*) is presented below. This alignment considers two incontinuous 3D substructures. Moreover, there is also incompatibility of chain id between compared 3D structures. Each substructure is described by id of its first residue [i.e., chain id + '\_' + residue serial number + '\_' + insertion code (if needed)] and length. To integrate many 3D substructures within a single alignment prepared for the particular RNA 3D structure(s) their descriptions are separated by '|'. Alignments prepared for the reference structure as well as the analyzed RNA 3D models combined into a single string are separated by ';'. Alignment prepared for the reference structure should be always included at the beginning of this string. 
```sh
solution.pdb:A_1,31|A_33,29;model.pdb:U_1,31|U_33,29
```
  
4. For all supported commands, a user can analyse both a single 3D model and many 3D models stored in the same directory at once.
  
### Important dependencies

RNAQUA uses the following external open source projects, namely:

- [Jersey][jersey] - framework for developing RESTful Web Services and their clients in Java that serves as a JAX-RS (JSR 311 & JSR 339) Reference Implementation,
- [Java Architecture for XML Binding][jaxb] - a library providing a fast and convenient way to bind XML schemas and Java representations to incorporate and process easily XML data,
- [Project Lombok][lombok] - a library allowing for compilation and building of a boilerplate-free code,
- [jarchivelib][jarchivelib] - an easy-to-use API layer on top of the [org.apache.commons.compress][org.apache.commons.compress].

RNAQUA is available in the [public repository][rnaqua] on GitHub.

### Requirements

To build the RNAQUA package one must have installed: 

- stable release of [Oracle JDK 7] [jdk] (recommended) or above, 
- stable release of [Apache Maven 3.2.3] [mvn] or above, 
- stable release of [Git] [git]. 

A used version of Java can be configured by setting the JAVA_HOME environment variable.

### Installation (common)

```sh
git clone https://github.com/mantczak/rnaqua.git rnaqua
cd rnaqua
```

### Build and tests (Windows)

```
build-and-tests.bat
```

### Build only (Windows)

```
build-only.bat
```

**_According to configuration of Linux/Mac machine (when maven3 package is installed, and 'No command mvn found') might be a need to add 'mvn3' symlink to 'mvn'._**

### Build and tests (Linux/Mac)

```sh
chmod u+x build-and-tests.sh
./build-and-tests.sh
```

### Build only (Linux/Mac)

```sh
chmod u+x build-only.sh
./build-only.sh
```

### Manual

The tool provides the following analysis modes: (1) a single 3D structure-based analysis, (2) analysis of 3D model(s) within the reference structure context.

- a single 3D structure-based analysis

```
 -c,--command <arg>                                  supported commands: PDB-VALIDATION, CLASH-SCORE, SEQUENCE, 
                                                     FRAGMENT, ORIGINAL-3D, RENUMERATED-3D
 -s,--single-model-file-path <arg>                   single model PDB file path
 -d,--multiple-models-directory-path <arg>           multiple PDB models directory path
 -o,--output-file-path <arg>                         output file path
 -a,--alignment <arg>                                (optional) single- or multi-models alignment
 -t,--base-pairs-identification-tool <arg>           (optional) base pairs identification tool, 
                                                     supported tools: RNAVIEW, MC-ANNOTATE 
                                                     [default=MC-ANNOTATE]
 -f,--consider-atoms-supported-by-RNA-Puzzles-only   (optional) consider atoms supported by RNA-Puzzles only
```

- analysis of 3D model(s) within the reference structure context
  
```
 -c,--command <arg>                                  supported commands: ROOT-MEAN-SQUARE-DEVIATION, 
                                                     ALL-INTERACTION-NETWORK-FIDELITY-SCORES-AT-ONCE,
                                                     INTERACTION-NETWORK-FIDELITY-WATSON-CRICK, 
                                                     INTERACTION-NETWORK-FIDELITY-NON-WATSON-CRICK,
                                                     INTERACTION-NETWORK-FIDELITY-STACKING, 
                                                     INTERACTION-NETWORK-FIDELITY-ALL, P-VALUE, 
                                                     DEFORMATION-INDEX, ALL-SCORES-AT-ONCE, SEQUENCE, FRAGMENT, 
                                                     ORIGINAL-3D, RENUMERATED-3D
 -r,--reference-structure-file-path <arg>            PDB file path of the reference structure
 -s,--single-model-file-path <arg>                   single model PDB file path
 -d,--multiple-models-directory-path <arg>           multiple PDB models directory path
 -o,--output-file-path <arg>                         output file path
 -a,--alignment <arg>                                (optional) single- or multiple-models alignment
 -t,--base-pairs-identification-tool <arg>           (optional) base pairs identification tool, 
                                                     supported tools: RNAVIEW, MC-ANNOTATE  
                                                     [default=MC-ANNOTATE]
 -f,--consider-atoms-supported-by-RNA-Puzzles-only   (optional) consider atoms supported by RNA-Puzzles only
 -h,--helices <arg>                                  (optional) helices description
 -l,--loops <arg>                                    (optional) loops description
 -w,--draw <arg>                                     (optional) draw description
 -m,--generate-matrix-data                           (optional) generate matrix data
 -g,--generate-svg-image                             (optional) generate svg image
```

### Tested configurations

- Linux Ubuntu 16.04 LTS x64, Oracle JDK 1.8.0_151 x64, Apache Maven 3.3.9.
- OS X El Capitan 10.11.3, Oracle JDK 1.7.0_80, Apache Maven 3.3.9.
- Windows 10 x64, Oracle JDK 1.7.0_80 i586, Apache Maven 3.2.3.

RNAQUA was tested on above configurations, but presumably it will work on other configurations too.

### Acknowledgements

We thank Prof. Eric Westhof and Zhichao Miao (RNA-Puzzles organizers) for sharing of ideas and discussions.

### Funding

The research was supported by the National Science Centre, Poland [grant No. 2012/06/A/ST6/00384].

License
----
Copyright (c) 2017 PUT Bioinformatics Group, licensed under [MIT license] [mit].

   [jersey]: https://jersey.github.io/
   [jaxb]: https://github.com/javaee/jaxb-v2
   [lombok]: https://projectlombok.org/
   [jarchivelib]: http://rauschig.org/jarchivelib/
   [org.apache.commons.compress]: http://commons.apache.org/proper/commons-compress/
   [jdk]: http://java.oracle.com/
   [mvn]: http://maven.apache.org/
   [git]: http://git-scm.com/
   [rnaqua]: https://github.com/mantczak/rnaqua.git
   [mit]: http://opensource.org/licenses/mit-license.php
   [molprobity]: http://scripts.iucr.org/cgi-bin/paper?S0907444909042073
   [rmsd]: http://scripts.iucr.org/cgi-bin/paper?S0567739476001873
   [inf]: https://dx.doi.org/10.1261%2Frna.1700409
   [pvalue]: https://dx.doi.org/10.1261%2Frna.1837410
   [rnassess]: https://doi.org/10.1093/nar/gkv557
