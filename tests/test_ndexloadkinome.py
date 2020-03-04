#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `ndexkinomeloader` package."""

import os
import tempfile
import shutil

import unittest
from ndexutil.config import NDExUtilConfig
from ndexkinomeloader import ndexloadkinome


class TestNdexkinomeloader(unittest.TestCase):
    """Tests for `ndexkinomeloader` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_parse_arguments(self):
        """Tests parse arguments"""
        res = ndexloadkinome._parse_arguments('hi', ['datadir'])

        self.assertEqual(res.profile, 'ndexkinomeloader')
        self.assertEqual(res.verbose, 0)
        self.assertEqual(res.logconf, None)
        self.assertEqual(res.conf, None)
        self.assertEqual(res.datadir, 'datadir')

        someargs = ['-vv','--conf', 'foo', '--logconf', 'hi',
                    '--profile', 'myprofy', 'datadir']
        res = ndexloadkinome._parse_arguments('hi', someargs)

        self.assertEqual(res.profile, 'myprofy')
        self.assertEqual(res.verbose, 2)
        self.assertEqual(res.logconf, 'hi')
        self.assertEqual(res.conf, 'foo')
        self.assertEqual(res.datadir, 'datadir')

    def test_setup_logging(self):
        """ Tests logging setup"""
        try:
            ndexloadkinome._setup_logging(None)
            self.fail('Expected AttributeError')
        except AttributeError:
            pass

        # args.logconf is None
        res = ndexloadkinome._parse_arguments('hi', ['datadir'])
        ndexloadkinome._setup_logging(res)

        # args.logconf set to a file
        try:
            temp_dir = tempfile.mkdtemp()

            logfile = os.path.join(temp_dir, 'log.conf')
            with open(logfile, 'w') as f:
                f.write("""[loggers]
keys=root

[handlers]
keys=stream_handler

[formatters]
keys=formatter

[logger_root]
level=DEBUG
handlers=stream_handler

[handler_stream_handler]
class=StreamHandler
level=DEBUG
formatter=formatter
args=(sys.stderr,)

[formatter_formatter]
format=%(asctime)s %(name)-12s %(levelname)-8s %(message)s""")

            res = ndexloadkinome._parse_arguments('hi', ['--logconf',
                                                         logfile,
                                                         'datadir'])
            ndexloadkinome._setup_logging(res)

        finally:
            shutil.rmtree(temp_dir)

    def test_main(self):
        """Tests main function"""

        temp_dir = tempfile.mkdtemp()
        # try where loading config is successful
        try:
            confile = os.path.join(temp_dir, 'some.conf')
            with open(confile, 'w') as f:
                f.write("""[hi]
                {user} = bob
                {pw} = smith
                {server} = dev.ndexbio.org""".format(user=NDExUtilConfig.USER,
                                                     pw=NDExUtilConfig.PASSWORD,
                                                     server=NDExUtilConfig.SERVER))
            res = ndexloadkinome.main(['myprog.py', '--conf',
                                       confile, '--profile',
                                       'hi', '--skipdownload',
                                       temp_dir])
            self.assertEqual(res, 2)
        finally:
            shutil.rmtree(temp_dir)
