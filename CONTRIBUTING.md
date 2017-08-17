## How to contribute

Contributions are much appreciated. The project follows the standard
[GitHub workflow](https://guides.github.com/introduction/flow/).

1. Fork the repository.
2. Create a new branch based on the `develop` branch (`git checkout -b your_branch develop`).
3. Do your modifications on that branch and commit your changes. Here's an
   example of a [good commit message](http://chris.beams.io/posts/git-commit/):
   ```
    [wip] Added Guo forcing scheme

    - Added Guo's getEqVelocityForcing method
    - Added Guo's getCollisionForcing method
    - Added Guo's getHydroVelocityForcing method
    - Added Guo option to the Create factory method in forcingScheme.h
    ```
4.  Push the changes to your fork (`git push origin your_branch`).
5.  Open a pull request against `develop`.

## Style guide

- Indent using 4 spaces.
- No trailing white spaces at the end of lines and no more than a single newline at the end of source file.
- Standard library `#include`s go first then a blank line followed by metaLBM's headers.
- Use your judgement but try to stick to the style of the surrounding code ([camelCase](https://en.wikipedia.org/wiki/Camel_case)).
