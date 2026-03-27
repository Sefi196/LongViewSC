# usage_logger.R
# Anonymous session logging to a shared SQLite database.
# No user-identifiable data is stored — only timestamps, app version,
# session duration, server RAM, and disk usage at session start.

USAGE_DB_PATH <- "/srv/shiny-server/usage_stats.db"

.init_usage_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), USAGE_DB_PATH)
  DBI::dbExecute(con, "PRAGMA journal_mode=WAL")
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS sessions (
      id            INTEGER PRIMARY KEY AUTOINCREMENT,
      app_version   TEXT    NOT NULL,
      started_at    TEXT    NOT NULL,
      ended_at      TEXT,
      duration_secs REAL,
      server_ram_mb REAL,
      disk_used_pct REAL
    )
  ")
  # Add disk_used_pct column if upgrading from older schema
  cols <- DBI::dbListFields(con, "sessions")
  if (!"disk_used_pct" %in% cols)
    DBI::dbExecute(con, "ALTER TABLE sessions ADD COLUMN disk_used_pct REAL")
  DBI::dbDisconnect(con)
}

.get_disk_pct <- function() {
  tryCatch({
    raw <- system("df /srv/shiny-server --output=pcent | tail -1", intern = TRUE)
    as.numeric(gsub("[^0-9]", "", trimws(raw)))
  }, error = function(e) NA_real_)
}

# Call at session start. Returns a session row id used by log_session_end().
log_session_start <- function(app_version) {
  tryCatch({
    .init_usage_db()
    disk_pct <- .get_disk_pct()
    con <- DBI::dbConnect(RSQLite::SQLite(), USAGE_DB_PATH)
    on.exit(DBI::dbDisconnect(con))
    DBI::dbExecute(con,
      "INSERT INTO sessions (app_version, started_at, disk_used_pct) VALUES (?, ?, ?)",
      list(app_version, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), disk_pct)
    )
    DBI::dbGetQuery(con, "SELECT last_insert_rowid() AS id")$id
  }, error = function(e) {
    warning("Usage log start failed: ", conditionMessage(e))
    NULL
  })
}

# Call in session$onSessionEnded().
log_session_end <- function(session_id, duration_secs, server_ram_mb) {
  if (is.null(session_id)) return(invisible(NULL))
  tryCatch({
    con <- DBI::dbConnect(RSQLite::SQLite(), USAGE_DB_PATH)
    on.exit(DBI::dbDisconnect(con))
    DBI::dbExecute(con,
      "UPDATE sessions SET ended_at=?, duration_secs=?, server_ram_mb=? WHERE id=?",
      list(
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        round(duration_secs, 1),
        round(server_ram_mb, 1),
        session_id
      )
    )
  }, error = function(e) {
    warning("Usage log end failed: ", conditionMessage(e))
  })
  invisible(NULL)
}
