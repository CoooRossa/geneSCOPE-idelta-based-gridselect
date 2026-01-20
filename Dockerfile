FROM rust:stable-slim AS build
WORKDIR /work
COPY . .
RUN cargo build --release

FROM debian:bookworm-slim
WORKDIR /app
COPY --from=build /work/target/release/idelta-gridselect /usr/local/bin/idelta-gridselect
COPY --from=build /work/config /app/config
ENTRYPOINT ["idelta-gridselect"]
