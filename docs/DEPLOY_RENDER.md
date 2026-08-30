# Deploy on Render (Standard)

Always-on hosting on the **Standard** compute plan (`1c-2g`, ~USD 25/month). This tier provides 2 GB RAM, which is required for whole-genome Biopython alignment in the analysis pipeline.

## One-time setup

1. Sign in at [render.com](https://render.com) (GitHub login works).
2. **New → Blueprint** (or **New → Web Service** if Blueprint is unavailable).
3. Connect repository: `coohzou/sars-cov-2-antibody-escape-predictor`.
4. Render reads `render.yaml` at the repo root:
   - **Plan:** `1c-2g` (Standard)
   - **Region:** Singapore (change in dashboard if you prefer Oregon/Frankfurt)
   - **Health check:** `/health`
5. Click **Apply** / **Create Web Service** and add a payment method when prompted.
6. Wait for the first build (typically 3–6 minutes).

Live URL:

`https://sars-cov-2-antibody-escape-predictor.onrender.com`

## Custom domain (optional)

1. Render dashboard → your service → **Settings → Custom Domains**.
2. Add your domain and follow the DNS CNAME instructions.
3. Render provisions HTTPS automatically.

## Verify deployment

```bash
curl https://YOUR-SERVICE.onrender.com/health
curl https://YOUR-SERVICE.onrender.com/ready
```

`/ready` should report `"predictor_ready": true` when models loaded correctly.

Upload a test FASTA (e.g. `data/prediction/gamma_complete.fasta`) from the web UI to confirm analysis completes without 502 errors.

## Updates

Push to `main`; Render redeploys automatically (`autoDeployTrigger: commit` in `render.yaml`).

## Paper / README

After the first successful deploy, set the live demo URL in your manuscript and in `README.md`.
