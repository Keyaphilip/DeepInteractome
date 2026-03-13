# Deployment Guide: Sharing DeepInteractome

The easiest way to share this application with your friend for free is to deploy it using a service called **Render**. Since we already committed the `Dockerfile`, Render can automatically build and host it for you!

Here are the exact steps to get it live:

## 1. Create a Render Account
1. Go to [Render (render.com)](https://render.com/) and sign up (you can use your GitHub account).

## 2. Connect Your Repository
1. In the Render Dashboard, click **New +** and select **Web Service**.
2. Select **"Build and deploy from a Git repository"**.
3. Connect your GitHub account if you haven't already.
4. Search for your repository `Keyaphilip/DeepInteractome` and click **Connect**.

## 3. Configure the Deployment
Fill in the deployment details as follows:
- **Name:** `deepinteractome` (or any name you like)
- **Region:** Pick whichever is closest to you.
- **Branch:** `main`
- **Root Directory:** (leave blank)
- **Environment:** `Docker` (Render should automatically detect the Dockerfile we made)
- **Instance Type:** Select the **Free** tier (this is perfect for sharing with your friend).

## 4. Advanced Settings (Important!)
Scroll down and click on **Advanced**. We need to tell Render what port to expect.
1. Click **Add Environment Variable**.
2. Add a variable with:
   - Key: `PORT`
   - Value: `8000`

## 5. Deploy
1. Scroll to the bottom and click **Create Web Service**.
2. Render will now pull your repository, build the Docker image, and deploy it. This usually takes a few minutes.

## 6. Share the Link!
Once the deployment finishes and the status turns to a green **Live**, you will see a URL near the top of the page (e.g., `https://deepinteractome-xxxx.onrender.com`).

Share this URL with your friend! 
*(Note: They will need to append `/ui` to the URL to access the web interface, e.g., `https://deepinteractome-xxxx.onrender.com/ui`)*
