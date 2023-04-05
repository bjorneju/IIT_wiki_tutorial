# Tutorial: Isolating Complexes and Unfolding $\Phi$-structures
---
This repository provides the supporting material for the "unfolding tutorial" within the [IIT wiki](https://centerforsleepandconsciousness.psychiatry.wisc.edu/intrinsic-ontology-wiki/). 
It is intended to help you understand what is happening "under the hood" in the process for isolating complexes and unfolding $\Phi$-structures within the integrated information theory (IIT). 
The tutorial is based on jupyter notebooks implementing functionality from [Pyphi](https://github.com/wmayner/pyphi/tree/feature/iit-4.0), but intentionally going into the details of the computational steps typically carried out behind the scenes. 
Keep in mind that the tutorial is best viewed as a part of the [IIT wiki](https://centerforsleepandconsciousness.psychiatry.wisc.edu/intrinsic-ontology-wiki/) as a whole. 
Also, if you are not interested in coding yourself, you can see screencasts of the tutorial content in the relevant pages of the wiki. 

For those of you interested in learning more, we strongly encourage you to head over to our wiki, read the academic papers (for example starting [here](https://arxiv.org/abs/2212.14787)), get acquainted with the latest version of pyphi (see documentation and examples for getting started [here](https://pyphi.readthedocs.io/en/latest/index.html)), browse other resources  online (for example [here](https://storage.googleapis.com/plos-corpus-prod/10.1371/journal.pcbi.1006343/2/pcbi.1006343.s001.pdf?X-Goog-Algorithm=GOOG4-RSA-SHA256&X-Goog-Credential=wombat-sa%40plos-prod.iam.gserviceaccount.com%2F20230308%2Fauto%2Fstorage%2Fgoog4_request&X-Goog-Date=20230308T091339Z&X-Goog-Expires=86400&X-Goog-SignedHeaders=host&X-Goog-Signature=1ddfbb8da0eb0b0cd31e7d3837e229d1abab6dcafabcd52032a10c4718d8100df6f3a07fa929323a2c884d01bcbd164c31b13bdf24557347e170721bbd216b5b7acd9aa19cf4d4447e1bec783f3ede8c2ee790c9c5ae9df1e71b4aff4c8198dc1d0db51650a97bcbcb7cbd0d7b074eac726d60c2e2862e98810e32caa7bada68f5484451decd96f5f55d15d0056cee819ffa85c2369a49dda03b3fa0fcf3997998da47a27621a2d850fda907d01c07b9eb647ec3691528a12bea02c05b03b74038335ac79ae8c8464ab7661c5e67a72f5ff503405b31f1b75dadd5d3593b4d0b031003b43968c42c6054c208f0e80f3c57b8a369bff3fe3e02e8ddd55e2b0923), [here](https://www.youtube.com/watch?v=0hex5katLGk&t=1630s) or [here](https://www.youtube.com/watch?v=FIZzxhJXJns)), or get in touch with us. 
However, be aware that a lot of the resources online are already outdated because the theory is still in development! 
And take a look at the our lab website ([the Center for Sleep and Consciousness](https://centerforsleepandconsciousness.psychiatry.wisc.edu/integrated-information-theory/)) for an update of the latest changes to the theory. 

---

To get started working with this tutorial on your own, please open your terminal and: 

1) make sure you have anaconda and git installed with a fresh base environment
2) in your terminal navigate to a folder where you would like to add the tutorial files
3) clone this repository with `git clone https://github.com/bjorneju/IIT_wiki_tutorial.git`
4) step into the folder using `cd IIT_wiki_tutorial`
5) create a new environment using `conda env create -f environment.yml` (this step takes some time)
6) activate your environment using `conda activate iit_tutorial`
7) open a jupyter lab session, and open the adress associated with the jupyter lab session (i.e. `http://localhost:8892/`)
8) and browse the notebooks!

Hope you enjoy the experience!
