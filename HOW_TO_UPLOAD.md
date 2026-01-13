# How to Upload EffOff to GitHub

Since the git and gh commands are not available in this environment, follow these manual steps to upload your package.

## 1. Create Repository on GitHub
1. Log in to [GitHub](https://github.com).
2. Click the **+** icon in the top-right and select **New repository**.
3. **Repository name**: `EffOff`
4. **Description**: Efficacy of Offset Forestry R Package
5. **Public/Private**: Public (required for easy install)
6. **Initialize this repository with**: [ ] Leave all unchecked (No README, No .gitignore, No License - we have these).
7. Click **Create repository**.

## 2. Push Local Code
Open your terminal (PowerShell, Command Prompt, or Git Bash) and run the following commands one by one:

```bash
# Navigate to the package directory
cd "A:\ColumbiaUniversity\AntiGravity_Workspaces\CAMPA\EffOff"

# Initialize Git
git init

# Add all files
git add .

# Commit files
git commit -m "Initial release of EffOff (v0.2.1)"

# Rename branch to main
git branch -M main

# Link to the GitHub repository you just created
git remote add origin https://github.com/vanubhav/EffOff.git

# Push to GitHub
git push -u origin main
```

**Note**: You may be prompted for your GitHub username and password. Use your Personal Access Token (PAT) as the password if 2FA is enabled.

## 3. Verify Installation
Once pushed, anyone can install the package using R:

```r
# install.packages("devtools")
devtools::install_github("vanubhav/EffOff")
```
