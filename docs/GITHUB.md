# Put this repository on GitHub

Repository title: **Black Swan Jeweled Math**

Repository name: `black-swan-jeweled-math`

Visibility: **Public**, as requested.

1. Create an empty repository at https://github.com/new?name=black-swan-jeweled-math . Leave the README, .gitignore, and license initialization options unchecked because this repository already contains its own files and initial commit.
2. Extract the prepared repository ZIP. Keep the complete directory structure, including its Git metadata.
3. Open a terminal in the extracted `black-swan-jeweled-math` directory and run:

```sh
git remote add origin https://github.com/ShadowAlzazel/black-swan-jeweled-math.git
git push -u origin main
```

The prepared repository already has an initial commit. You do not need to run git init or make a duplicate initial commit.

If you create a repository with a different owner or name, substitute the URL GitHub shows you. Authenticate with your own GitHub account when Git asks; the archive contains no credentials.

## Future changes

```sh
git add .
git commit -m "Describe the new experiment"
git push
```

## Optional web gallery

The explorers are standalone HTML and can be opened from the web directory without publishing. To publish the journey page and the labs:

1. Open the repository's Settings → Pages.
2. Under Build and deployment, choose Source: Deploy from a branch; Branch: `main`, folder `/ (root)`; Save.
3. After the first deployment finishes (a minute or two), the journey is at https://shadowalzazel.github.io/black-swan-jeweled-math/ and the lab gallery at https://shadowalzazel.github.io/black-swan-jeweled-math/web/ .

The root `index.html` is the front page, `.nojekyll` serves every file as-is, and link previews use `images/social_card.png`. Pushing alone does not publish anything until Pages is turned on.
