FROM registry.ccb-gitlab.cs.uni-saarland.de/pascal/aviator/aviator_production_django

FROM nginx:latest
COPY compose/production/nginx/nginx.conf /etc/nginx/nginx.conf
COPY --from=0 /app/staticfiles /static
