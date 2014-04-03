/*
 * Run the following commands as the MySQL superuser to create the
 * UCSC user and database.
 */
create user 'ucsc'@'localhost' identified by 'ucsc';
create database ucsc;
grant all on ucsc.* to 'ucsc'@'localhost';
grant file on *.* to 'ucsc'@'localhost';
flush privileges;
