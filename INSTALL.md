# <img src="https://www.docker.com/wp-content/uploads/2022/03/Moby-logo.png" alt="Docker Logo" height="75"/>
<img src="https://www.nicepng.com/png/full/506-5068487_open-install-png-logo.png" alt="Install" height="20"/> First Install [Docker](https://docs.docker.com/get-docker/)

**On macOS and Linux**

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src="https://icon-library.com/images/play-icon-white-png/play-icon-white-png-26.jpg" alt="Run" height="20"/> `docker run --rm -v $PWD:$PWD -w $PWD joan/3dna python -m 3dna data/plasmid_8k.fasta`


**On Windows**
   1. open a PowerShell and run:

   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src="https://icon-library.com/images/play-icon-white-png/play-icon-white-png-26.jpg" alt="Run" height="20"/> `wsl --update`
   
   2. Restart your computer

   3. Run Docker Desktop

   4. In a PowerShell:

   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src="https://icon-library.com/images/play-icon-white-png/play-icon-white-png-26.jpg" alt="Run" height="20"/> `docker run --rm -v <absolute_path_with_/_instead_of_\>:/3dna joan/3dna bash -c 'cd /3dna ; python -m 3dna data/plasmid_8k.fasta'`

<br/>
<br/>

# <img src="https://upload.wikimedia.org/wikipedia/commons/e/ea/Conda_logo.svg" alt="Conda Logo" height="50"/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src="https://www.nicepng.com/png/full/506-5068487_open-install-png-logo.png" alt="Install" height="20"/> Install first [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src="https://www.nicepng.com/png/full/506-5068487_open-install-png-logo.png" alt="Install" height="20"/> `conda env create -f environment.yml`

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Then,

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src="https://icon-library.com/images/play-icon-white-png/play-icon-white-png-26.jpg" alt="Run" height="20"/> `conda run -n 3dna python -m 3dna data/plasmid_8k.fasta`

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;or

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src="https://icon-library.com/images/play-icon-white-png/play-icon-white-png-26.jpg" alt="Run" height="20"/> `conda activate 3dna`

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<img src="https://icon-library.com/images/play-icon-white-png/play-icon-white-png-26.jpg" alt="Run" height="20"/> `python -m dna data/plasmid_8k.fasta`
