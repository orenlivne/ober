/*
 * Run the following commands as the MySQL 'impute' user to delete and recreate the
 * 'impute' database schema.
 */
create user 'impute'@'localhost' identified by 'impute';
create database impute;
grant all on impute.* to 'impute'@'localhost';
flush privileges;
