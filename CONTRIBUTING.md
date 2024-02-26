# Contributing to volmdlr

First off, thanks for taking the time to contribute! 

All types of contributions are encouraged and valued. The community looks forward to your contributions. 

> And if you like the project, but just don't have time to contribute, that's fine. There are other easy ways to support the project and show your appreciation, which we would also be very happy about:
> - Use volmdlr in your projects!
> - Star the project
> - Talk about it
> - Refer this project in your project's readme
> - Mention the project at local meetups and tell your friends/colleagues


## I want to contribute to the code

### Ask before coding

If you are not sure about your contribution you can ask the question on the issues or MP the project leader.

### Process

We have a process based on three branches:

- dev that receive the features and impacting changes pull requests
- testing that is regularly a freeze of testing. Release candidate (RC) can be made on testing
- master that receive the code from testing before the release. Release are made on master

If you want to contribute, you need first to fork this project on your repository, following the GitHub directive. See:
https://git-scm.com/book/en/v2/GitHub-Contributing-to-a-Project

Once you have forked the project you can create a thematic branch where you will add your contribution and then open a
pull request

### CI

Your contribution may be refused by our CI. In this case:
- lowering the coverage: add some tests
- introducing pylint errors
- triggering pydocstyle errors
- non-respect of PEP8
- no CHANGELOG.md edit

Correct the CI errors to get approval from volmdlr team!

## I Have a Question

Before you ask a question, it is best to search for existing [Issues](/issues) that might help you. In case you have found a suitable issue and still need clarification, you can write your question in this issue. It is also advisable to search the internet for answers first.

If you then still feel the need to ask a question and need clarification, we recommend the following:

- Open an [Issue](/issues/new).
- Provide as much context as you can about what you're running into.
- Provide project and platform versions , depending on what seems relevant.
- Put the label "Question" to your issue
