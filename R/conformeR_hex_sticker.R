library(hexSticker)
library(ggplot2)
library(showtext)

font_add_google("Roboto", "roboto")
showtext_auto()
showtext_opts(dpi = 300)   # must match the dpi used in sticker()

set.seed(42)

# ---- Shape helpers ------------------------------------------
circle <- function(x0, y0, r, n = 100) {
  th <- seq(0, 2 * pi, length.out = n)
  data.frame(x = x0 + r * cos(th), y = y0 + r * sin(th))
}
ellipse <- function(x0, y0, rx, ry, angle = 0, from = 0, to = 2 * pi, n = 100) {
  th <- seq(from, to, length.out = n)
  x <- rx * cos(th); y <- ry * sin(th)
  data.frame(x = x0 + x * cos(angle) - y * sin(angle),
             y = y0 + x * sin(angle) + y * cos(angle))
}
rect <- function(x1, x2, y1, y2) data.frame(x = c(x1, x2, x2, x1), y = c(y1, y1, y2, y2))
poly <- function(x, y) data.frame(x = x, y = y)
ink  <- "#1b1b1b"
shape <- function(df, fill, col = ink, lw = 0.25) {
  geom_polygon(data = df, aes(x, y), fill = fill, colour = col,
               linewidth = lw, inherit.aes = FALSE)
}

# ---- Muted, heraldic palette ----------------------------------
pal     <- c("#b5523b", "#3d7f78", "#c9a24a", "#6b7f5a")  # clusters
fur     <- "#8f8a84"
fur_lt  <- "#e8e4dc"
steel   <- "#a7b0b9"
steel_d <- "#5c6670"
gold    <- "#b8973f"
navy    <- "#1f2f4a"
wine    <- "#7a1e24"
amber   <- "#d98e04"

# ---- Background clusters ------------------------------------
centers <- data.frame(x = c(-0.66, 0.66, -0.74, 0.78),
                      y = c( 0.55, 0.78, -0.22, -0.55))

clusters <- do.call(rbind, lapply(seq_len(nrow(centers)), function(i) {
  do.call(rbind, lapply(1:12, function(j) {
    d <- circle(centers$x[i] + rnorm(1, 0, 0.07),
                centers$y[i] + rnorm(1, 0, 0.07),
                runif(1, 0.02, 0.035), n = 40)
    d$id <- paste(i, j); d$cluster <- factor(i)
    d
  }))
}))

hulls <- do.call(rbind, lapply(seq_len(nrow(centers)), function(i) {
  d <- circle(centers$x[i], centers$y[i], 0.22)
  d$cluster <- factor(i)
  d
}))

# ---- Ringed tail ------------------------------------------------
tail_bands <- function(N = 20, w0 = 0.06, w1 = 0.04) {
  t  <- seq(0, 1, length.out = N + 1)
  x  <- 0.22 + 0.55 * sin(t * pi * 0.6)
  y  <- -0.60 + 1.10 * t
  dx <- c(diff(x), diff(x)[N]); dy <- c(diff(y), diff(y)[N])
  len <- sqrt(dx^2 + dy^2); nx <- -dy / len; ny <- dx / len
  w  <- seq(w0, w1, length.out = N + 1)
  lx <- x + nx * w; ly <- y + ny * w
  rx <- x - nx * w; ry <- y - ny * w
  bands <- do.call(rbind, lapply(1:N, function(k) {
    data.frame(x = c(lx[k], lx[k + 1], rx[k + 1], rx[k]),
               y = c(ly[k], ly[k + 1], ry[k + 1], ry[k]),
               band = k,
               col = ifelse(((k - 1) %/% 2) %% 2 == 0, "light", "dark"))
  }))
  list(bands = bands, tip = circle(x[N + 1], y[N + 1], w1))
}
tl <- tail_bands()

# ---- Illustration -----------------------------------------------
p <- ggplot() +
  # clusters
  geom_path(data = hulls, aes(x, y, group = cluster, colour = cluster),
            linetype = "dashed", linewidth = 0.2, alpha = 0.7) +
  geom_polygon(data = clusters, aes(x, y, group = id, fill = cluster),
               alpha = 0.6, colour = NA) +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +

  # tail
  geom_polygon(data = subset(tl$bands, col == "light"), aes(x, y, group = band),
               fill = fur_lt, colour = ink, linewidth = 0.15) +
  geom_polygon(data = subset(tl$bands, col == "dark"), aes(x, y, group = band),
               fill = ink, colour = ink, linewidth = 0.15) +
  shape(tl$tip, ink) +

  # sabatons, body, pauldrons
  shape(ellipse(-0.13, -0.76, 0.10, 0.045), steel_d) +
  shape(ellipse( 0.13, -0.76, 0.10, 0.045), steel_d) +
  shape(ellipse(0, -0.38, 0.30, 0.35), steel) +
  shape(ellipse(-0.24, -0.10, 0.10, 0.07), steel) +
  shape(ellipse( 0.24, -0.10, 0.10, 0.07), steel) +

  # surcoat with gold trim and belt
  shape(poly(c(-0.20, 0.20, 0.26, -0.26), c(-0.08, -0.08, -0.62, -0.62)),
        wine, col = gold, lw = 0.4) +
  shape(rect(-0.23, 0.23, -0.44, -0.40), steel_d) +
  shape(rect(-0.025, 0.025, -0.445, -0.395), gold) +

  # shield: chevron + cluster charges
  shape(poly(c(-0.55, -0.11, -0.11, -0.15, -0.22, -0.33, -0.44, -0.51, -0.55),
             c(-0.15, -0.15, -0.40, -0.52, -0.63, -0.72, -0.63, -0.52, -0.40)),
        navy, col = gold, lw = 0.6) +
  shape(poly(c(-0.53, -0.33, -0.13, -0.13, -0.33, -0.53),
             c(-0.42, -0.26, -0.42, -0.50, -0.34, -0.50)), gold, col = NA) +
  shape(circle(-0.40, -0.58, 0.026), pal[1], col = NA) +
  shape(circle(-0.33, -0.58, 0.026), pal[2], col = NA) +
  shape(circle(-0.26, -0.58, 0.026), pal[3], col = NA) +

  # ears
  shape(poly(c(-0.28, -0.19, -0.35), c(0.28, 0.33, 0.45)), fur) +
  shape(poly(c( 0.28,  0.19,  0.35), c(0.28, 0.33, 0.45)), fur) +

  # head, face, snout
  shape(ellipse(0, 0.21, 0.24, 0.20), fur) +
  shape(ellipse(0, 0.17, 0.16, 0.10), fur_lt, col = NA) +
  shape(poly(c(-0.07, 0.07, 0), c(0.15, 0.15, 0.06)), fur_lt, col = NA) +
  # almond eyes in sloping dark patches
  shape(ellipse(-0.08, 0.20, 0.058, 0.040, angle = -0.5), ink, col = NA) +
  shape(ellipse( 0.08, 0.20, 0.058, 0.040, angle =  0.5), ink, col = NA) +
  shape(ellipse(-0.08, 0.20, 0.030, 0.019, angle = -0.5), amber, col = NA) +
  shape(ellipse( 0.08, 0.20, 0.030, 0.019, angle =  0.5), amber, col = NA) +
  shape(ellipse(-0.08, 0.20, 0.007, 0.016), ink, col = NA) +
  shape(ellipse( 0.08, 0.20, 0.007, 0.016), ink, col = NA) +
  # nose, placed higher on the snout
  shape(ellipse(0, 0.105, 0.028, 0.018), ink, col = NA) +

  # helmet: dome, rim, nasal guard, rivets
  shape(ellipse(0, 0.30, 0.26, 0.23, from = 0, to = pi), steel) +
  shape(rect(-0.005, 0.005, 0.30, 0.53), steel_d, col = NA) +          # ridge
  shape(rect(-0.28, 0.28, 0.27, 0.32), steel_d) +
  shape(rect(-0.013, 0.013, 0.20, 0.30), steel_d, lw = 0.15) +
  shape(circle(-0.20, 0.295, 0.009), steel, col = NA) +
  shape(circle(-0.10, 0.295, 0.009), steel, col = NA) +
  shape(circle( 0.10, 0.295, 0.009), steel, col = NA) +
  shape(circle( 0.20, 0.295, 0.009), steel, col = NA) +

  # sword + armoured arm
  shape(rect(0.408, 0.432, -0.20, 0.42), "#d5dbe0") +
  shape(poly(c(0.408, 0.432, 0.42), c(0.42, 0.42, 0.48)), "#d5dbe0") +
  shape(rect(0.418, 0.422, -0.18, 0.38), steel_d, col = NA) +          # fuller
  shape(rect(0.35, 0.49, -0.23, -0.205), gold) +
  shape(rect(0.408, 0.432, -0.35, -0.23), "#4e342e") +
  shape(circle(0.42, -0.37, 0.022), gold) +
  shape(ellipse(0.34, -0.30, 0.08, 0.045), steel) +
  shape(circle(0.42, -0.29, 0.038), steel_d) +                          # gauntlet

  coord_fixed(xlim = c(-1, 1), ylim = c(-0.9, 0.95), expand = FALSE) +
  theme_void() +
  theme_transparent() +
  theme(legend.position = "none")

# ---- Hex sticker --------------------------------------------------
sticker(
  p,
  package  = "conformeR",
  p_size   = 7,
  p_family = "roboto",
  p_color  = gold,
  p_y      = 1.48,
  s_x      = 1,
  s_y      = 0.74,
  s_width  = 1.15,
  s_height = 1.05,
  h_fill   = "#1b1f2a",
  h_color  = gold,
  h_size   = 1.2,
  filename = "conformeR_hex.png",
  dpi      = 300
)

usethis::use_logo("conformeR_hex.png")
