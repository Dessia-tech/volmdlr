# Contributing to volmdlr

Thank you for considering contributing to volmdlr! Your contributions are highly appreciated and valued by the community.

> If you find the project valuable, but you are unable to contribute directly, there are other ways you can show your support:
> - Use volmdlr in your projects!
> - Star the project on GitHub.
> - Spread the word about the project.
> - Include a reference to the project in your project's README file.
> - Mention the project during local meetups and share it with your friends and colleagues.


## I want to contribute to the code

### Ask before coding

If you're uncertain about your contribution, feel free to reach out by opening an issue or sending a message to the project leader.

### Process

We have a process based on three branches:

- dev that receive the features and impacting changes pull requests
- testing that is regularly a freeze of testing. Release candidate (RC) can be made on testing
- master that receive the code from testing before the release. Release are made on master

To contribute code:

1. Fork the project to your GitHub repository following the instructions [here.](https://git-scm.com/book/en/v2/GitHub-Contributing-to-a-Project)
2. Create a thematic branch where you'll add your contribution.
3. Open a pull request once you've added your contribution to your branch, ensuring that you select the 'dev' branch as the target branch in volmdlr repository.

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
