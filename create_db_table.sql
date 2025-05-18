CREATE TABLE `loggerInfo` (
  `idx` bigint NOT NULL AUTO_INCREMENT,
  `loggerID` varchar(11) DEFAULT NULL,
  `latitude` double DEFAULT NULL,
  `longitude` double DEFAULT NULL,
  `loggerPosition` text,
  `available` double DEFAULT NULL,
  `bathymetry` double DEFAULT NULL,
  `year` bigint DEFAULT NULL,
  `site` varchar(11) DEFAULT NULL,
  `notedDepth` float DEFAULT NULL,
  KEY `ix_loggerInfo_index` (`idx`)
) ENGINE=InnoDB AUTO_INCREMENT=844 DEFAULT CHARSET=latin1;


CREATE TABLE `loggerData_2021` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `Time` datetime DEFAULT NULL,
  `DO` float DEFAULT NULL,
  `Temp` float DEFAULT NULL,
  `logger` varchar(11) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Time` (`Time`),
  KEY `logger` (`logger`)
) ENGINE=InnoDB AUTO_INCREMENT=2333766 DEFAULT CHARSET=utf8mb3;



CREATE TABLE `loggerData_2022` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `Time` datetime DEFAULT NULL,
  `DO` float DEFAULT NULL,
  `Temp` float DEFAULT NULL,
  `logger` varchar(11) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Time` (`Time`),
  KEY `logger` (`logger`)
) ENGINE=InnoDB AUTO_INCREMENT=2333766 DEFAULT CHARSET=utf8mb3;


CREATE TABLE `loggerData_2023` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `Time` datetime DEFAULT NULL,
  `DO` float DEFAULT NULL,
  `Temp` float DEFAULT NULL,
  `logger` varchar(11) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `Time` (`Time`),
  KEY `logger` (`logger`)
) ENGINE=InnoDB AUTO_INCREMENT=2333766 DEFAULT CHARSET=utf8mb3;