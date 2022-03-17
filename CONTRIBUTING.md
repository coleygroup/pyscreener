# How to contribute

We welcome contributions from external contributors, and this document
describes how to merge code changes into `pyscreener`. 

## Getting Started

1. Make sure you have a [GitHub account](https://github.com/signup/free).
1. [Fork](https://help.github.com/articles/fork-a-repo/) this repository on GitHub.
1. On your local machine, [clone](https://help.github.com/articles/cloning-a-repository/) your fork of the repository.

## Making Changes

1. Install the development version of `pyscreener`: `pip install -e .[dev]`
1. Install the the `pre-commit` hooks: `pre-commit install`
1. Add some changes to your local fork.  It's usually a [good idea](http://blog.jasonmeridth.com/posts/do-not-issue-pull-requests-from-your-master-branch/) to make changes on a [branch](https://help.github.com/articles/creating-and-deleting-branches-within-your-repository/) with the branch name relating to the feature you are going to add.
1. If you're providing a new feature, you **must** add test cases and documentation.
1. Complete the checklist below
1. When you are ready for others to examine and comment on your new feature, navigate to your fork of `pyscreener` on GitHub and open a [pull request](https://help.github.com/articles/using-pull-requests/) (PR). Note that after you launch a PR from one of your fork's branches, all subsequent commits to that branch will be added to the open pull request automatically.  Each commit added to the PR will be validated for mergability, compilation and test suite compliance; the results of these tests will be visible on the PR page.
1. When you believe the PR has reached its final draft, check the "Ready to go" box below to let the `pyscreener` devs know that the changes are complete. The code will not be considered for merging until this box is checked, the continuous integration returns checkmarks, and multiple core developers give "Approved" reviews.

## Checklist

- [ ] formatted with black: `black pyscreener`
- [ ] linted with flake8: `flake8`
- [ ] unit tests added (if necessary)?
- [ ] all unit tests pass (if any core code was changed)?

## Ready to merge?

- [ ] yes!

# Additional Resources

* [General GitHub documentation](https://help.github.com/)
* [PR best practices](http://codeinthehole.com/writing/pull-requests-and-other-good-practices-for-teams-using-github/)
* [A guide to contributing to software packages](http://www.contribution-guide.org)
* [Thinkful PR example](http://www.thinkful.com/learn/github-pull-request-tutorial/#Time-to-Submit-Your-First-PR)
