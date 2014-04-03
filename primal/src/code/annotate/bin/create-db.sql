-- Hutterites CGI variant annotation database
-- Create database and account
--
-- Run with:
-- sudo -i
-- su - postgres
-- cat create-db.sql | psql

CREATE database hutt if not exists;
CREATE USER hutt WITH PASSWORD 'hutt' if not exists;
GRANT ALL PRIVILEGES ON DATABASE 'hutt' to hutt;
