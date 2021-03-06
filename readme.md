NSAPH Project Template
================
Ben Sabath
04/10/2019

This directory is an premade structure that can be directly cloned from
github to serve as a starting structure for new projects done as part of
the team. This document will provide basic instructions on how to use
this template to get started on your project.

# Downloading this repository

You can copy this repository directly from github to a directory of your
choice using the following
command:

``` bash
git clone https://github.com/NSAPH/project_template.git <NEW_PROJECT_NAME>
```

You will likely need to set up an rsa key on your code.harvard.edu
account to use the remote git features. You can find information on how
to set that up
[here](https://help.github.com/enterprise/2.12/user/articles/adding-a-new-ssh-key-to-your-github-account).
You will need at least one key set up for each computer you are working
on.

TODO:

Confirm the location (github vs code.harvard.edu)

# Setting up Git

After initially cloning the repository, the downloaded repository will
be connected to the template repository. In order to set up your own
project, you will first need to disconnect the repository from the
project template remote repository. Use the following command:

    git remote rm origin

You will then need to set up your own repository for the project. First,
create an empty repository (make sure to not include the readme.md
option that github suggests) for your project, either on github or on
code.harvard.edu. We recommend using code.harvard.edu, as it will make
all materials only available within Harvard rather than fully open to
the public. If you would like to make your repository public at a later
date, copying it from code.harvard.edu to github is fairly
straightforward.

Next, run the following lines to connect the template to the new
repository:

Connecting to a code.harvard.edu
    repository:

    git remote add origin git@code.harvard.edu:<username>/<repository_name>.git
    git push -u origin master

Connecting to a github repository:

    git remote add origin https://github.com/<username>/<project>.git
    git push -u origin master

You can then use git as you normally would.

# Other initial setup

The readmes (including this file) will pertain to the project template
not to your project. You should edit them (and remove extra files) so
that they describe what the intent of your research project is.

# The .gitignore file

The following is the default .gitignore file included in the directory.

``` r
cat(scan(".gitignore", what="character"), sep = "\n")
```

    ## data/*.csv
    ## data/*.rds
    ## data/*.RDS
    ## results/*.rds
    ## results/*.RDS

The .gitignore file tells git which files it should ignore. Files listed
in the .gitignore will not have their changes tracked and will not be
uploaded to the remote server. This is useful as we want to avoid
tracking large files (such as data files, which are also frequently
sensitive and shouldn’t be stored on less secure systems like github).
If there are other files you would like to avoid tracking, edit the
.gitignore file to add their name.

# Directories

  - `code`
      - The code directory should be used to store all code used as part
        of the project. If there are a large number of code files
        additional directories within the `code` directory are
        recommended, especially if there are multiple workflows. There
        should also be a description of the order code should be run in,
        or some other description of the contents and a method of
        indicating the workflow.
  - `data`
      - This directory should contain all raw data and processed data
        used by a given project. All csv and rds files within this
        directory are assumed to be large and by default are not tracked
        by git. If there are a large number of data files additional
        internal structure is recommended. The `.gitignore` file will
        need to be updated if more structure is used. One common
        paradign is to have a `data/raw_data` directory storing unedited
        files as receinved from the source and a `data/analysis_data`
        folder containing the assembled and cleaned data ready for use
        with models.
  - `figures`
      - This directory should contain all of the figures generated
        through the course of reaserch. The readme file in this
        directory should list the files, have a brief discription of the
        figures, and list the file that creates the figure. This
        directory should be tracked by git.
  - `reports`
      - This directory should contain Rmarkdown files used to summarize
        and describe reserach processes and data features. The HTML or
        PDF outputs of the markdown files should also be stored in this
        directory. This directory should be tracked by git.
  - `results`
      - This directory should be used to store objects such as models,
        tables, or other similar products of analysis. Whether or not
        files in this directory should be tracked by github is a
        question of their size and sensitivity and likely varies by
        project. By default, .RDS files in this directory will not be
        tracked, but all other files will.
