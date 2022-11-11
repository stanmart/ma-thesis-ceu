# ma-thesis-ceu

This repo contains the codes for my MA thesis at the Central European University. You can find the actual paper [here](https://www.etd.ceu.edu/2016/stancsics_martin.pdf).


## Compilation
Disclaimer: This project was created before I knew how to do properly reproducible research. This repo (or, more precisely, its [Bitbucket version](https://bitbucket.org/stanmart/thesis/src/master/) was just meant to be a quick way to share the codes without too many instructions. Anyways, if you would like to reproduce my results, here are some pointers:


  * The folder `huggett` contains the codes used in the thesis. The Jupyter notebook in that folder is meant to be used as a script you run to generate the results, while the `.py` files in the same folder contain functions imported by the notebook.
  * You will need the packages listed in the `requirements.txt` file. It's anyone's guess though which version of those packages and Python itself was used at the time. My best guess is those included in some Anaconda 4.* distribution with Python 3 from around 2016 spring. I would probably start with something like the following:
    * `python==3.5.2`
    * `notebook=4.2.1`
    * `numpy==1.11.1`
    * `scipy==0.17.1`
    * `matplotlib==1.5.1`
    * `seaborn==0.7.0`
  * The folder `pre-alpha` contains some very incomplete Haskell code for calculating the model in a more complex setting. It is just very early attempt, does not work, and was not used for the thesis.
